
#ifndef REV_ANALYSIS_H
#define REV_ANALYSIS_H
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
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <fstream>

#include "common/ScaLBL.h"
#include "common/Communication.h"
#include "common/MPI.h"
#include "analysis/Minkowski.h"
#include "ProfilerApp.h"
#include "models/MRTModel.h"



class REVfunc {
public:
    double average_pore_size;
    double rev_x_poro, rev_x_perm, rev_x_tort, rev_x_surf;
    double size_step;
    int det_samples, stat_samples;
    double gamma;
    void PoreSize(ScaLBL_MRTModel& model);
    void DetRevAnalysis(ScaLBL_MRTModel& model, string filename);
    void StatRevAnalysis(ScaLBL_MRTModel& model, string filename);
};

#endif // REV_ANALYSIS_H
