#STATUS:
#!bin/bash

main() {
    get_info_for_experiment
    determine_set_of_controls_based_on_sample_data
    grab_and_output_to_to_testData
}
