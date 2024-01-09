#!/bin/bash

echo "Not stable yet. Nothing happened."
exit 1

EMA_PATH="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
TEMP_DIR="$(mktemp -d)"

DO_NOT_DELETE_FILE_LIST_IN_REPORTS_DIR=("style.css" "model.html")
DO_NOT_DELETE_FILE_LIST_IN_GRAPHICS_DIR=()

mv $EMA_PATH/share/graphics/* $TEMP_DIR
for i in "${DO_NOT_DELETE_FILE_LIST_IN_GRAPHICS_DIR[@]}"; do
    mv $TEMP_DIR/$i $EMA_PATH/share/graphics
done

mv $EMA_PATH/share/reports/* $TEMP_DIR
for i in "${DO_NOT_DELETE_FILE_LIST_IN_REPORTS_DIR[@]}"; do
    mv $TEMP_DIR/$i $EMA_PATH/share/reports
done

echo "All useless files have been moved to : $TEMP_DIR. The temporary directory will be deleted when you will log out of the current session."
# revert operation ?
