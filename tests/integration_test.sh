TROPOMI_DIR=/data/uptrop/nobackup/tropomi/Data/
PANDORA_DIR=/data/uptrop/nobackup/pandora/
GC_DIR=/data/uptrop/Projects/DEFRA-NH3/GC/
OUT_DIR=./integration_outputs

echo Testing tropomi_ut_no2...
python ../uptrop/tropomi_ut_no2.py --trop_dir $TROPOMI_DIR --out_dir $OUT_DIR \
--season djf --grid_res 4x5 --cloud_product fresco --cloud_threshold 07

echo Testing ut_no2_gc_test.py
python ../uptrop/ut_no2_gc_test.py --gc_dir $GC_DIR --out_dir $OUT_DIR

echo Testing compare_tropomi_pandora.py
python ../uptrop/compare_tropomi_pandora.py --trop_dir $TROPOMI_DIR --pan_dir $PANDORA_DIR --out_dir $OUT_DIR \
--start_date 2019-06-01 --end_date 2019-06-30 --apply_bias_correction True

echo Testing fresco_cl_err.py
python ../uptrop/fresco_cld_err.py --s5p_data_dir $TROPOMI_DIR --output_dir $OUT_DIR
