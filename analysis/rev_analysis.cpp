/*
  Copyright 2013--2018 James E. McClure, Virginia Polytechnic & State University
  Copyright Equnior ASA

  This file is part of the Open Porous Media project (OPM).
  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*
 * Multi-relaxation time LBM Model
 */
#include "analysis/rev_analysis.h"
#include "models/MRTModel.h"
#include "analysis/distance.h"
#include "common/ReadMicroCT.h"
#include <cmath>
#include <algorithm>
#include <random>


struct det_convergence_metrics {
    double RE;
    double CC;
};

struct stat_convergence_metrics {
    double CC, mean, max, min;
};

struct VoxelDist {
    int x, y, z;
    double dist;
};

det_convergence_metrics det_compute_convergence(const std::vector<double>& data, int iter_step, int det_samples) {
    det_convergence_metrics m = {-1, -1};
    
    if (iter_step < 1) 
        return m;

    double mean = 0.0;
    double std_dev = 0.0;
    int valid_count = 0;
    int start_idx = (iter_step >= det_samples - 1) ? (iter_step + 1 - det_samples) : 0;
    for (int l = start_idx; l <= iter_step; l++) {
        if (data[l] != -1) {
            mean += data[l];
            valid_count++;
        }
    }

    if (valid_count == 0) return m;
    
    mean /= valid_count;

    for (int l = start_idx; l <= iter_step; l++) {
        if (data[l] != -1) {
            std_dev += (data[l] - mean) * (data[l] - mean);
        }
    }
    std_dev = sqrt(std_dev / valid_count);

    if (mean != 0) {
        m.CC = std_dev / mean;
    }

    if (data[iter_step] != -1 && data[iter_step - 1] != -1) {
        double sum_prev_curr = data[iter_step] + data[iter_step - 1];
        if (sum_prev_curr != 0) {
            m.RE = 2.0 * fabs((data[iter_step] - data[iter_step - 1]) / sum_prev_curr);
        }
    }

    return m;
}

stat_convergence_metrics stat_compute_convergence(const std::vector<double>& data, int samples) {
    stat_convergence_metrics res = {-1.0, 0, -1, 1e100};
    int valid_count = 0;

    for (int l = 0; l < samples; l++) {
        if (data[l] != -1) {
            res.mean += data[l];
            if (data[l] > res.max) res.max = data[l];
            if (data[l] < res.min) res.min = data[l];
            valid_count++;
        }
    }

    if (valid_count == 0) return res;

    res.mean /= valid_count;

    double std_dev = 0.0;
    for (int l = 0; l < samples; l++) {
        if (data[l] != -1) {
            std_dev += (data[l] - res.mean) * (data[l] - res.mean);
        }
    }
    std_dev = std::sqrt(std_dev / valid_count);
    
    if (res.mean != 0) res.CC = std_dev / res.mean;

    return res;
}

void REVfunc::PoreSize(ScaLBL_MRTModel &MRT) {
    int Nx = MRT.Nx, Ny = MRT.Ny, Nz = MRT.Nz;
    double h = MRT.Dm->voxel_length;

    DoubleArray distance_updated;
    distance_updated.resize(Nx, Ny, Nz);
    
    std::vector<VoxelDist> pore_voxels;

    for (int k = 0; k < Nz; k++) {
        for (int j = 0; j < Ny; j++) {
            for (int i = 0; i < Nx; i++) {
                double val = MRT.Distance(i, j, k);
                distance_updated(i, j, k) = val;
                
                if (i > 0 && i < Nx - 1 && j > 0 && j < Ny - 1 && k > 0 && k < Nz - 1) {
                    if (val > 0) {
                        pore_voxels.push_back({i, j, k, val});
                    }
                }
            }
        }
    }

    std::sort(pore_voxels.begin(), pore_voxels.end(), [](const VoxelDist& a, const VoxelDist& b) {
        return a.dist > b.dist;
    });

    for (const auto& voxel : pore_voxels) {
        int i = voxel.x;
        int j = voxel.y;
        int k = voxel.z;
        double current_max = voxel.dist;

        int start_x = std::max(1, static_cast<int>(i - current_max - 1));
        int end_x   = std::min(Nx - 1, static_cast<int>(i + current_max + 2));
        int start_y = std::max(1, static_cast<int>(j - current_max - 1));
        int end_y   = std::min(Ny - 1, static_cast<int>(j + current_max + 2));
        int start_z = std::max(1, static_cast<int>(k - current_max - 1));
        int end_z   = std::min(Nz - 1, static_cast<int>(k + current_max + 2));

        double r_sq = current_max * current_max;

        for (int n = start_z; n < end_z; n++) {
            for (int m = start_y; m < end_y; m++) {
                for (int l = start_x; l < end_x; l++) {
                    double dist_sq = (i - l) * (i - l) + (j - m) * (j - m) + (k - n) * (k - n);
                    
                    if (dist_sq <= r_sq) {
                        double value2 = distance_updated(l, m, n);
                        if (value2 > 0 && value2 < current_max) {
                            distance_updated(l, m, n) = current_max;
                        }
                    }
                }
            }
        }
    }

    double sum = 0.0;
    int total_pore_voxels = 0;
    
    std::map<double, int> psd_histogram; 

    for (int k = 1; k < Nz - 1; k++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int i = 1; i < Nx - 1; i++) {
                double val = distance_updated(i, j, k);
                if (val > 0) {
                    sum += val;
                    total_pore_voxels++;
                    psd_histogram[val]++;
                }
            }
        }
    }

    average_pore_size = (total_pore_voxels > 0) ? (2.0 * sum / total_pore_voxels) : 0.0;

    std::ofstream outfile("pore_size_distribution.csv");
    if (outfile.is_open()) {
        outfile << "radius_voxels diameter_um voxel_count volume_fraction\n";
        
        for (const auto& pair : psd_histogram) {
            double radius_vx = pair.first;
            int count = pair.second;
            
            double diameter_um = 2.0 * radius_vx * h;
            double volume_fraction = static_cast<double>(count) / total_pore_voxels;
            
            outfile << radius_vx << " "
                    << diameter_um << " " 
                    << count << " " 
                    << volume_fraction << "\n";
        }
        
        outfile << "\nAverage_Pore_Size_um " << average_pore_size * h << "\n";
        outfile.close();
    }

    std::cout << "\n\n\nAverage Pore Size: " << average_pore_size * h << " micrometers\n\n\n" << std::endl;
}

void REVfunc::DetRevAnalysis(ScaLBL_MRTModel &MRT, string filename) {

    auto db = std::make_shared<Database>(filename);
    auto mrt_db = db->getDatabase("MRT");

    bool run_rev = false;
    if(mrt_db->keyExists("REV"))
        run_rev = mrt_db->getScalar<bool>("REV");

    if(run_rev == false)
        return;


    size_step = 0.1;
    det_samples = 10;
    gamma = 0.1;


    auto rev_db = db->getDatabase("REV");
    if(rev_db->keyExists("SizeStep"))
        size_step = rev_db->getScalar<double>("SizeStep");
    if(rev_db->keyExists("DetSamples"))
        det_samples = rev_db->getScalar<int>("DetSamples");
    if(rev_db->keyExists("gamma"))
        gamma = rev_db->getScalar<double>("gamma");
        
    PoreSize(MRT);

    double current_size_x_double = 2.0 * average_pore_size;
    int x_size = static_cast<int>(std::ceil(current_size_x_double)); //could be a problem if the geometry is in the yz plane.


    int Nx = MRT.Nx, Ny = MRT.Ny, Nz = MRT.Nz;
    int iter_step = 0; 
    double h = MRT.Dm->voxel_length;
    double mu = (MRT.tau - 0.5) / 3.f;
    
    rev_x_poro = -1;
    rev_x_perm = -1;
    rev_x_tort = -1;
    rev_x_surf = -1;

    int vector_size_x = 1000;
    std::vector<double> sub_x_size(vector_size_x);
    std::vector<double> perm(vector_size_x);
    std::vector<double> poro(vector_size_x);
    std::vector<double> tort(vector_size_x);
    std::vector<double> surf(vector_size_x);

    std::ofstream log_file("det_rev_analysis.csv", std::ios::app);
    if (log_file.tellp() == 0) {
    log_file << "iter_step x_size_um x_size_vx "
             << "porosity rev_poro_um rev_poro_vx RE_poro CC_poro "
             << "permeability rev_perm_um rev_perm_vx RE_perm CC_perm "
             << "tortuosity rev_tort_um rev_tort_vx RE_tort CC_tort "
             << "surface rev_surf_um rev_surf_vx RE_surf CC_surf\n";
    }


    while (x_size < Nx - 2) {

        //(still have to test for 2d geometries)
        int y_size = (Ny - 2) * x_size / (Nx - 2);
        int z_size = (Nz - 2) * x_size / (Nx - 2);

        int start_x = (Nx - 2 - x_size) / 2 + 1;
        int start_y = (Ny - 2 - y_size) / 2 + 1;
        int start_z = (Nz - 2 - z_size) / 2 + 1;

        int end_x = (Nx - 2 + x_size) / 2 + 1;
        int end_y = (Ny - 2 + y_size) / 2 + 1;
        int end_z = (Nz - 2 + z_size) / 2 + 1;

        start_x = std::max(1, start_x);
        end_x = std::min(Nx-2, end_x);
        start_y = std::max(1, start_y);
        end_y = std::min(Ny-2, end_y);
        start_z = std::max(1, start_z);
        end_z = std::min(Nz-2, end_z);

        double current_surf = 0;
        double count = 0.0;
        double vax = 0.0, vay = 0.0, vaz = 0.0;
        double v_tort=0.0;

        for (int k = start_z; k <= end_z; k++) { 
            for (int j = start_y; j <= end_y; j++) {
                for (int i = start_x; i <= end_x; i++) {
                    if (MRT.Distance(i, j, k) > 0) {
                        vax += MRT.Velocity_x(i, j, k);
                        vay += MRT.Velocity_y(i, j, k);
                        vaz += MRT.Velocity_z(i, j, k);
                        v_tort += sqrt(MRT.Velocity_x(i, j, k)*MRT.Velocity_x(i, j, k) + MRT.Velocity_y(i, j, k)*MRT.Velocity_y(i, j, k) + MRT.Velocity_z(i, j, k)*MRT.Velocity_z(i, j, k));
                        count += 1.0;
                    }
                    if ((MRT.Distance(i, j, k) * MRT.Distance(i, j, k + 1) < 0) && (k != end_z))
                        current_surf++;

                    if ((MRT.Distance(i, j, k) * MRT.Distance(i, j + 1, k) < 0) && (j != end_y))
                        current_surf++;

                    if ((MRT.Distance(i, j, k) * MRT.Distance(i + 1, j, k) < 0) && (i != end_x))
                        current_surf++;
                }
            }
        }

        if (count == 0)
        {
            vax = 0.0;
            vay = 0.0;
            vaz = 0.0;
        } else {
        vax /= count;
        vay /= count;
        vaz /= count;}

        double force_mag = sqrt(MRT.Fx * MRT.Fx + MRT.Fy * MRT.Fy + MRT.Fz * MRT.Fz);
        double dir_x = MRT.Fx / force_mag;
        double dir_y = MRT.Fy / force_mag;
        double dir_z = MRT.Fz / force_mag;
        if (force_mag == 0.0) {
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
            force_mag = 1.0;
        }

        current_surf = current_surf / ((end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1));
        double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);
        double current_poro = count / ((end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1));
        double current_perm = 1013 * h * h * mu * current_poro * flow_rate / force_mag;
        double current_tort = (vaz != 0) ? (v_tort / (vaz * count) - 1.0) : -1;

        sub_x_size[iter_step] = x_size;
        poro[iter_step] = current_poro;
        perm[iter_step] = current_perm;
        tort[iter_step] = current_tort;
        surf[iter_step] = current_surf;

        //Porosity deterministic approach
        det_convergence_metrics poro_m = det_compute_convergence(poro, iter_step, det_samples);
        if (poro_m.RE >= 0 && poro_m.RE < gamma && poro_m.CC >= 0 && poro_m.CC < gamma && rev_x_poro == -1) {
            rev_x_poro = x_size;
        }
        log_file << iter_step << " " << x_size * h << " " << x_size << " "
                 << poro[iter_step] << " " << rev_x_poro*h << " " << rev_x_poro << " " << poro_m.RE << " " << poro_m.CC << " ";

        //Permeability deterministic approach
        det_convergence_metrics perm_m = det_compute_convergence(perm, iter_step, det_samples);
        if (perm_m.RE >= 0 && perm_m.RE < gamma && perm_m.CC >= 0 && perm_m.CC < gamma && rev_x_perm == -1) {
            rev_x_perm = x_size;
        }
        log_file << perm[iter_step] << " " << rev_x_perm*h << " " << rev_x_perm << " " << perm_m.RE << " " << perm_m.CC << " ";

        //Tortuosity deterministic approach
        det_convergence_metrics tort_m = det_compute_convergence(tort, iter_step, det_samples);
        if (tort_m.RE >= 0 && tort_m.RE < gamma && tort_m.CC >= 0 && tort_m.CC < gamma && rev_x_tort == -1) {
            rev_x_tort = x_size;
        }
        log_file << tort[iter_step]+1 << " " << rev_x_tort*h << " " << rev_x_tort << " " << tort_m.RE << " " << tort_m.CC << " ";

        //Surface Deterministic approach
        det_convergence_metrics surf_m = det_compute_convergence(surf, iter_step, det_samples);
        if (surf_m.RE >= 0 && surf_m.RE < gamma && surf_m.CC >= 0 && surf_m.CC < gamma && rev_x_surf == -1) {
            rev_x_surf = x_size;
        }
        log_file << surf[iter_step] << " " << rev_x_surf*h << " " << rev_x_surf << " " << surf_m.RE << " " << surf_m.CC << "\n";

        double min_step = std::max(1.0, average_pore_size);
        double proportional_step = current_size_x_double * size_step;
        current_size_x_double += std::max(min_step, proportional_step);
        x_size = static_cast<int>(std::ceil(current_size_x_double)); 

        iter_step++;
    }

        //whole domain data

        double current_surf = 0;
        double vax = 0.0, vay = 0.0, vaz = 0.0;
        double count = 0.0;
        double v_tort = 0.0;

        for (int k = 1; k <= Nz-2; k++) { 
            for (int j = 1; j <= Ny-2; j++) {
                for (int i = 1; i <= Nx-2; i++) {
                    if (MRT.Distance(i, j, k) > 0) {
                        vax += MRT.Velocity_x(i, j, k);
                        vay += MRT.Velocity_y(i, j, k);
                        vaz += MRT.Velocity_z(i, j, k);
                        v_tort += sqrt(MRT.Velocity_x(i, j, k)*MRT.Velocity_x(i, j, k) + MRT.Velocity_y(i, j, k)*MRT.Velocity_y(i, j, k) + MRT.Velocity_z(i, j, k)*MRT.Velocity_z(i, j, k));
                        count += 1.0;
                    }
                    if ((MRT.Distance(i, j, k) * MRT.Distance(i, j, k + 1) < 0) && (k != Nz-2))
                        current_surf++;

                    if ((MRT.Distance(i, j, k) * MRT.Distance(i, j + 1, k) < 0) && (j != Ny-2))
                        current_surf++;

                    if ((MRT.Distance(i, j, k) * MRT.Distance(i + 1, j, k) < 0) && (i != Nx-2))
                        current_surf++;
                }
            }
        }

        if (count == 0)
        {
            vax = 0.0;
            vay = 0.0;
            vaz = 0.0;
        } else {
        vax /= count;
        vay /= count;
        vaz /= count;}

        double force_mag = sqrt(MRT.Fx * MRT.Fx + MRT.Fy * MRT.Fy + MRT.Fz * MRT.Fz);
        double dir_x = MRT.Fx / force_mag;
        double dir_y = MRT.Fy / force_mag;
        double dir_z = MRT.Fz / force_mag;
        if (force_mag == 0.0) {
            dir_x = 0.0;
            dir_y = 0.0;
            dir_z = 1.0;
            force_mag = 1.0;
        }

        double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);
        double current_poro = count / ((Nx-2) * (Ny-2) * (Nz-2));
        double current_perm = 1013 * h * h * mu * current_poro * flow_rate / force_mag;
        double current_tort = (vaz != 0) ? (v_tort / (vaz * count) - 1.0) : -1;
        current_surf = current_surf / ((Nx-2) * (Ny-2) * (Nz-2));

        log_file << "-1 " << (Nx-2)*h << " " << (Nx-2) << " " 
                 << current_poro << " " << rev_x_poro*h << " " << rev_x_poro << " 0 0 " 
                 << current_perm << " " << rev_x_perm*h << " " << rev_x_perm << " 0 0 " 
                 << current_tort << " " << rev_x_tort*h << " " << rev_x_tort << " 0 0 " 
                 << current_surf << " " << rev_x_surf*h << " " << rev_x_surf << " 0 0\n";

            log_file.close();
            printf("\n\nEnded Deterministic Analysis \n\n");
}

void REVfunc::StatRevAnalysis(ScaLBL_MRTModel &MRT, string filename) {

    std::mt19937 rng(0);

    auto db = std::make_shared<Database>(filename);
    auto mrt_db = db->getDatabase("MRT");
    
    if (!mrt_db->keyExists("REV") || !mrt_db->getScalar<bool>("REV"))
        return;

    stat_samples = 50;
    auto rev_db = db->getDatabase("REV");
    if (rev_db->keyExists("StatSamples"))
        stat_samples = rev_db->getScalar<int>("StatSamples");

    int Nx = MRT.Nx, Ny = MRT.Ny, Nz = MRT.Nz;
    
    double h = MRT.Dm->voxel_length;
    double mu = (MRT.tau - 0.5) / 3.0f;

    std::vector<double> stat_poro(stat_samples);
    std::vector<double> stat_perm(stat_samples);
    std::vector<double> stat_tort(stat_samples);
    std::vector<double> stat_surf(stat_samples);

    std::ofstream log_file("stat_rev_analysis.csv", std::ios::app);
    if (log_file.tellp() == 0) {
        log_file << "x_size_um x_size_vx stat_poro stat_poro_mean stat_poro_max stat_poro_min "
                 << "stat_perm stat_perm_mean stat_perm_max stat_perm_min "
                 << "stat_tort stat_tort_mean stat_tort_max stat_tort_min "
                 << "stat_surf stat_surf_mean stat_surf_max stat_surf_min\n";
    }

    double current_size_x_double = 1e100;
    if (rev_x_poro > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_poro);
    if (rev_x_perm > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_perm);
    if (rev_x_tort > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_tort);
    if (rev_x_surf > 0) current_size_x_double = std::min(current_size_x_double, (double)rev_x_surf);
    if (current_size_x_double == 1e100) {
        current_size_x_double = 2.0 * average_pore_size;
    }
    if (current_size_x_double == -1) current_size_x_double = 2.0 * average_pore_size;
    int x_size = static_cast<int>(std::ceil(current_size_x_double));

    while (x_size > 0 && x_size < Nx - 2) {
        
        int sub_y = std::floor(x_size * (Ny - 2) / (Nx - 2)) + 1;
        int sub_z = std::floor(x_size * (Nz - 2) / (Nx - 2)) + 1;

        int range_x = (Nx - 2) - x_size + 1;
        int range_y = (Ny - 2) - sub_y + 1;
        int range_z = (Nz - 2) - sub_z + 1;

        if (range_x < 1) range_x = 1;
        if (range_y < 1) range_y = 1;
        if (range_z < 1) range_z = 1;

        std::uniform_int_distribution<int> dist_x(1, range_x);
        std::uniform_int_distribution<int> dist_y(1, range_y);
        std::uniform_int_distribution<int> dist_z(1, range_z);

        for (int l = 0; l < stat_samples; l++) {
            int start_x = dist_x(rng);
            int start_y = dist_y(rng);
            int start_z = dist_z(rng);

            int end_x = std::min(Nx - 2, start_x + x_size - 1);
            int end_y = std::min(Ny - 2, start_y + sub_y - 1);
            int end_z = std::min(Nz - 2, start_z + sub_z - 1);

            double vax = 0.0, vay = 0.0, vaz = 0.0, count = 0.0, v_tort = 0.0, current_surf = 0.0;

            for (int k = start_z; k <= end_z; k++) {
                for (int j = start_y; j <= end_y; j++) {
                    for (int i = start_x; i <= end_x; i++) {
                        double dist_val = MRT.Distance(i, j, k);
                        if (dist_val > 0) {
                            vax += MRT.Velocity_x(i, j, k);
                            vay += MRT.Velocity_y(i, j, k);
                            vaz += MRT.Velocity_z(i, j, k);
                            v_tort += std::sqrt(std::pow(MRT.Velocity_x(i, j, k), 2) + 
                                                std::pow(MRT.Velocity_y(i, j, k), 2) + 
                                                std::pow(MRT.Velocity_z(i, j, k), 2));
                            count += 1.0;
                        }
                        
                        if (k != end_z && (dist_val * MRT.Distance(i, j, k + 1) < 0)) current_surf++;
                        if (j != end_y && (dist_val * MRT.Distance(i, j + 1, k) < 0)) current_surf++;
                        if (i != end_x && (dist_val * MRT.Distance(i + 1, j, k) < 0)) current_surf++;
                    }
                }
            }

            if (count > 0) {
                vax /= count; vay /= count; vaz /= count;
            }

            double force_mag = std::sqrt(MRT.Fx * MRT.Fx + MRT.Fy * MRT.Fy + MRT.Fz * MRT.Fz);
            double dir_x = (force_mag > 0) ? MRT.Fx / force_mag : 0.0;
            double dir_y = (force_mag > 0) ? MRT.Fy / force_mag : 0.0;
            double dir_z = (force_mag > 0) ? MRT.Fz / force_mag : 1.0;
            if (force_mag == 0.0) force_mag = 1.0;

            double vol = (end_x - start_x + 1) * (end_y - start_y + 1) * (end_z - start_z + 1);
            current_surf /= vol;
            
            double flow_rate = (vax * dir_x + vay * dir_y + vaz * dir_z);
            double current_poro = count / vol;
            double current_perm = 1013.0 * h * h * mu * current_poro * flow_rate / force_mag;
            double current_tort = (vaz != 0) ? (v_tort / (vaz * count) - 1.0) : -1;

            stat_poro[l] = current_poro;
            stat_perm[l] = current_perm;
            stat_tort[l] = current_tort;
            stat_surf[l] = current_surf;
        }

        stat_convergence_metrics p_stats = stat_compute_convergence(stat_poro, stat_samples); 
        stat_convergence_metrics k_stats = stat_compute_convergence(stat_perm, stat_samples); 
        stat_convergence_metrics t_stats = stat_compute_convergence(stat_tort, stat_samples); 
        stat_convergence_metrics s_stats = stat_compute_convergence(stat_surf, stat_samples); 
        log_file << x_size * h << " " << x_size << " "
                << p_stats.CC << " " << p_stats.mean << " " << p_stats.max << " " << p_stats.min << " "
                << k_stats.CC << " " << k_stats.mean << " " << k_stats.max << " " << k_stats.min << " "
                << t_stats.CC << " " << t_stats.mean + 1 << " " << t_stats.max + 1 << " " << t_stats.min + 1 << " "
                << s_stats.CC << " " << s_stats.mean / h << " " << s_stats.max / h << " " << s_stats.min / h << "\n";

        double min_step = std::max(1.0, average_pore_size);
        double proportional_step = current_size_x_double * size_step;
        current_size_x_double += std::max(min_step, proportional_step);
        x_size = static_cast<int>(std::ceil(current_size_x_double));
    }
    
    log_file.close();
    printf("\n\nEnded Statistical Analysis\n\n");
}