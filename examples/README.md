### How to make your own automatic test on TITAN:

- Start a new empty folder in `TITAN/examples`. This folder will be used in `.gitlab-ci.yml` in the `TEST_FOLDER` variable.
- Copy/create only the necessary files: `input`, `basis`, elemental files and maybe `kbands` and `inputmag`.
- If there is no self-consistency for this system:
  * Run a calculation with `-> itype=1`
  * `cp -r results results_correct`
  * `sed -n '/Self-consistent ground state/,/\(.*iterations \)/p' output/* | head -n -1 | cut -c-77 > results_correct/scResult`
  * Add a job for the self-consistency in the `.gitlab-ci.yml` following the pattern there
  * Push the changes and check if the pipeline will work (otherwise, work on the fix)
- For a new calculation (`itype>1`)
  * run the calculation with the desired itype
  * `cp -r results results_correct`
  * Add a job in `.gitlab-ci.yml` following the pattern, with the corresponding self-consistency as dependency

*Notes:*
- Have in mind that good tests are finished in ~<1min. They don’t need to have physical meaning or accuracy, so the number of points can be used to speed it up.
- It is recommended to run a local test (a second and maybe even third run) after copying the files to `results_correct`, and do the local comparison with:
```for file in `find results/ -type f` ; do echo "Comparig file $(basename ${file})" ; python ../../scripts/compare_output.py ${file} ${file/results/results_correct} ; done```
- If there are temporary files that should not be compared, use the EXCLUDE variable in `.gitlab-ci.yml`
- The default tolerance (`TOL=1.d-7`) is preferable. If it needs to be increased, maybe something is wrong in the calculation.
