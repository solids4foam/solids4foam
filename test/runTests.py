##---------------------------------------------------------------------------##
##    Description
##        This file is part of the solids4foam test suite.
##
##        For every sub-directory in the "testCases" directory which contains a
##        "testExpectedValues.json" file, copy the contents of the equivalent
##        tutorial case from "../tutorials" to this directory. Then execute the
##        Allrun script and perform the values checks in the
##        "testExpectedValues.json" file against values in the generated
##        "log.solids4Foam" file. The "testExpectedValues.json" file is given in
##        following format:
##
##        {
##           "log_pattern": "Max sigmaEq \\(von Mises stress\\) = ([\\d\\.]+)",
##           "comparisons": [
##             {
##               "type": "equality",
##               "value": 1.62914
##             }
##           ]
##        }
##
##        where the "log_pattern" is the pattern to be extracted from
##        "log.solids4Foam" and "comparisons" can be "equality", "range",
##        "greater_than" or "less_than". See "exampleTestExpectedValues.json"
##        for detailed usage.
##
##        Execute this script from the current directory with:
##            > python runTests.py
##
##        You can clean the testCases directory with:
##            > ./Allclean
##
##    Author
##        Philip Cardiff, UCD.
##
##---------------------------------------------------------------------------##
import subprocess
import os
import re
import sys
import json
import shutil

def copy_directory_contents(src, dst):
    try:
        # Ensure the destination directory exists
        if not os.path.exists(dst):
            os.makedirs(dst)

        # Iterate over the contents of the source directory
        for item in os.listdir(src):
            s = os.path.join(src, item)
            d = os.path.join(dst, item)
            
            if os.path.isdir(s):
                # Recursively copy subdirectory contents
                shutil.copytree(s, d)
            else:
                # Copy file
                shutil.copy2(s, d)
        
        print(f"Contents of {src} copied to {dst}")
    except shutil.Error as e:
        print(f"Error: {e}")
    except OSError as e:
        print(f"OS error: {e}")

def run_allrun(script_path):
    script_dir = os.path.dirname(script_path)
    script_name = os.path.basename(script_path)
    result = subprocess.run(['bash', script_name], cwd=script_dir, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return result.returncode, result.stdout, result.stderr

def parse_log(log_file, pattern):
    values = []
    regex = re.compile(pattern)
    with open(log_file, 'r') as file:
        for line in file:
            match = regex.search(line)
            if match:
                values.append(float(match.group(1)))
    return values

def load_expected_values(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

def compare_values(generated_values, comparison):
    comparison_type = comparison["type"]
    expected_value = comparison["value"]
    
    for value in generated_values:
        if comparison_type == "equality":
            if value != expected_value:
                return False, f"Value mismatch: expected {expected_value}, got {value}"
        elif comparison_type == "range":
            min_val, max_val = expected_value
            if not (min_val <= value <= max_val):
                return False, f"Value out of range: expected {min_val} <= value <= {max_val}, got {value}"
        elif comparison_type == "greater_than":
            if value <= expected_value:
                return False, f"Value not greater than expected: expected > {expected_value}, got {value}"
        elif comparison_type == "less_than":
            if value >= expected_value:
                return False, f"Value not less than expected: expected < {expected_value}, got {value}"
    return True, ""

def run_tests(test_dir):
    for root, dirs, files in os.walk(test_dir):
        for dir in dirs:
            test_case_dir = os.path.join(root, dir)
            allrun_script = os.path.join(test_case_dir, 'Allrun')
            log_file = os.path.join(test_case_dir, 'log.solids4foam')
            expected_values_file = os.path.join(test_case_dir, 'testExpectedValues.json')
            
            #if os.path.exists(allrun_script) and os.path.exists(expected_values_file):
            if os.path.exists(expected_values_file):
                print(f"\nRunning test: {test_case_dir}")

                # Copy the contents of the appropriate tutorial case
                tut_rel_path = test_case_dir[len("./testCases/"):]
                tut_dir = os.path.join('../tutorials', tut_rel_path)
                copy_directory_contents(tut_dir, test_case_dir)

                returncode, stdout, stderr = run_allrun(allrun_script)
                if returncode != 0:
                    print(f"Test {test_case_dir} failed: Allrun script returned code {returncode}")
                    print("stderr:")
                    print(stderr.decode())
                    print("stdout:")
                    print(stdout.decode())
                    continue
                
                expected_values = load_expected_values(expected_values_file)
                log_pattern = expected_values["log_pattern"]
                comparisons = expected_values["comparisons"]
                
                generated_values = parse_log(log_file, log_pattern)
                
                for comparison in comparisons:
                    success, message = compare_values(generated_values, comparison)
                    if not success:
                        print(f"Test {test_case_dir} FAILED: {message}")
                        break
                else:
                    print(f"Test {test_case_dir} PASSED.")
            #else:
            #    print(f"Skipping test {test_case_dir}: Missing Allrun or testExpectedValues.json")

def main():
    test_dir = './testCases'
    run_tests(test_dir)

if __name__ == "__main__":
    main()
