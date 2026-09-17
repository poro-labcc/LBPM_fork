import os
import glob
import numpy as np
import pyvista as pv

def load_vtk_fields(filepath):
    """
    Reads a VTK file (legacy, XML, or Parallel) using PyVista 
    and extracts the fq_XX fields as NumPy arrays.
    """
    try:
        # PyVista natively handles .pvti, .vti, and legacy .vtk
        mesh = pv.read(filepath)
    except Exception as e:
        print(f"    [!] Failed to read {filepath}: {e}")
        return None

    fields = {}
    for q in range(19):
        field_name = f"fq_{q:02d}"
        
        # Check cell data first, then point data
        if field_name in mesh.cell_data:
            fields[field_name] = np.array(mesh.cell_data[field_name])
        elif field_name in mesh.point_data:
            fields[field_name] = np.array(mesh.point_data[field_name])
        else:
            print(f"    [!] Warning: Field '{field_name}' not found in {filepath}")

    return fields

def detailed_vtk_analysis(base_filepath, other_filepath):
    print(f"    [*] Loading and parsing VTK data with PyVista...")
    base_fields = load_vtk_fields(base_filepath)
    other_fields = load_vtk_fields(other_filepath)

    if not base_fields or not other_fields:
        print("    [FAIL] Could not load fields from one or both files.")
        return

    max_prints = 10
    total_diffs = 0

    for q in range(19):
        field_name = f"fq_{q:02d}"
        if field_name not in base_fields or field_name not in other_fields:
            continue

        base_arr = base_fields[field_name]
        other_arr = other_fields[field_name]

        if base_arr.shape != other_arr.shape:
            print(f"    [FAIL] Size mismatch in {field_name}: Base is {base_arr.shape}, Other is {other_arr.shape}")
            continue

        # Find indices where values differ by more than a tiny tolerance
        # (Avoids false positives from MPI floating-point truncation)
        mismatches = np.where(np.abs(base_arr - other_arr) > 1e-14)[0]
        
        if len(mismatches) > 0:
            total_diffs += len(mismatches)
            print(f"    [!] {field_name} has {len(mismatches)} mismatched voxels.")
            
            # Print the first few mismatches for this specific field
            for idx in mismatches[:max_prints]:
                diff = abs(base_arr[idx] - other_arr[idx])
                print(f"      -> {field_name}, Global Voxel {idx:8d} | Base: {base_arr[idx]: .8e} | Other: {other_arr[idx]: .8e} | Delta: {diff:.2e}")
            
            if len(mismatches) > max_prints:
                print(f"      -> ... and {len(mismatches) - max_prints} more in {field_name}\n")

    if total_diffs == 0:
        print("    [OK]   All 19 fields match perfectly!")
    else:
        print(f"    [!] Total differences found across all channels: {total_diffs}")


def compare_vtk_debug(folders):
    if len(folders) < 2:
        print("Error: Provide at least two folders to compare.")
        return

    base_folder = folders[0]
    
    # Point directly to the summary.pvti
    search_pattern_pvti = os.path.join(base_folder, "summary.pvti")
    
    base_files = glob.glob(search_pattern_pvti)

    if not base_files:
        print(f"No VTK debug files found matching '{search_pattern_pvti}'")
        return

    all_match = True

    print(f"Found {len(base_files)} VTK debug file(s) in base folder '{base_folder}'.\n")

    for base_filepath in sorted(base_files):
        filename = os.path.basename(base_filepath)
        print(f"--- Checking {filename} ---")

        for other_folder in folders[1:]:
            # Removed the hardcoded "vis" since the folders array already contains it
            other_filepath = os.path.join(other_folder, filename)

            if not os.path.exists(other_filepath):
                print(f"  [FAIL] Missing in '{other_folder}'")
                all_match = False
                continue

            # Perform the mathematical comparison of the VTK arrays
            print(f"  Comparing '{base_folder}' vs '{other_folder}'...")
            detailed_vtk_analysis(base_filepath, other_filepath)

    print("\n==================================================")
    if all_match:
        print("FINISHED COMPARISON: SUCCESS")
    else:
        print("FINISHED COMPARISON: FAILURES DETECTED")
    print("==================================================")

# =================================================================
# DEFINE YOUR FOLDERS HERE
# =================================================================
target_folders = [
    "./DEBUG_PARALLEL_INIT/1_1_1/vis000/",
    "./DEBUG_PARALLEL_INIT/1_1_2/vis000/",
    "./DEBUG_PARALLEL_INIT/1_2_1/vis000/",
    "./DEBUG_PARALLEL_INIT/2_1_1/vis000/",
    "./DEBUG_PARALLEL_INIT/2_2_2/vis000/",
]

compare_vtk_debug(target_folders)