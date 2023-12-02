#!/bin/sh

set -u

#### This script will run Next flow script models_mutations.groovy

./nextflow run ../models_mutations.groovy --root_dir ./ --target_batch_limit 2
