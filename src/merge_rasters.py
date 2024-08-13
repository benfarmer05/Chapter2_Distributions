import os
import rasterio
from rasterio.merge import merge
import numpy as np

def read_rasters(file_list):
    datasets = []
    for file in file_list:
        src = rasterio.open(file)
        datasets.append(src)
    return datasets

def merge_rasters(datasets):
    # Merge the rasters giving priority to the highest resolution
    merged_data, merged_transform = merge(datasets, method='first')
    merged_meta = datasets[0].meta.copy()
    merged_meta.update({
        "driver": "GTiff",
        "height": merged_data.shape[1],
        "width": merged_data.shape[2],
        "transform": merged_transform
    })
    return merged_data, merged_meta

def save_raster(output_path, data, meta):
    with rasterio.open(output_path, 'w', **meta) as dest:
        dest.write(data)

def main():
    # List of TIFF files to be merged
    tiff_files = [
        'path/to/your/first.tiff',
        'path/to/your/second.tiff',
        'path/to/your/third.tiff'
        # Add more TIFF file paths as needed
    ]

    # Read the rasters
    datasets = read_rasters(tiff_files)

    # Merge the rasters
    merged_data, merged_meta = merge_rasters(datasets)

    # Save the merged raster
    output_path = 'path/to/output/merged.tiff'
    save_raster(output_path, merged_data, merged_meta)

    # Close all datasets
    for dataset in datasets:
        dataset.close()

if __name__ == "__main__":
    main()