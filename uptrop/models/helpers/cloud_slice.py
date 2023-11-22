from .cloud_slice_no2 import cloud_slice_no2
from .cloud_slice_o3 import cloud_slice_o3

CLOUD_SLICE_ERROR_ENUM = {
    1: "too_few_points",
    2: "low_cloud_height_range",
    3: "low_cloud_height_std",
    4: "large_error",
    5: "sig_diff_from_zero",
    6: "outlier",
    7: "non_uni_strat"
}


# cloud_slice_no2 and cloud_slice_o3 are allowed to evolve independently and take different input arguments
# But we may not need this, if the processors are built differently using cloud_slice_no2 or cloud_slice_o3
def cloud_slice(species, *args, **kwargs):    
    if species.lower() == "no2":
        return cloud_slice_no2(*args, **kwargs)
    elif species.lower() == "o3":
        return cloud_slice_o3(*args, **kwargs)
    else:
        valid_species = ["no2", "o3"]
        valid_species_upper = [s.upper() for s in valid_species]
        raise ValueError(
            f"Invalid species '{species}'. Supported values: {', '.join(valid_species + valid_species_upper)}")