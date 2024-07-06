#!/bin/bash

create_new_experiment() {
    local date_of_creation=$(date +%Y%m%d)
    local filename="${date_of_creation}_experiment.md"
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

}
create_new_experiment
