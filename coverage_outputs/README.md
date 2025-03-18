# Generating Code coverage
At the moment this is slightly stone-henge, due to the fact that CUDA is not supported for code coverage, and you can't automate a GPU test for free. For now, users must run the unit tests locally (using `ctest`, see [ReadTheDocs](https://woden.readthedocs.io/en/latest/testing/cmake_testing.html) for instructions). They can then run `source create_cov_reports.sh`, which will grab the outputs of the `C` tests from `ctest` and covert them into an appropriate formats, as well as running the `python` tests using `coverage`, which also generates appropriate outputs. To get `coverage`, do something like:

```bash
pip install coverage
```

Once those outputs have been created, you can then run `source send_reports_to_codecov.sh` to update the [codecov](https://about.codecov.io/) hosted coverage report. You will need an environment variable `WODEN_CODECOV_TOKEN` to be able to do this. Read the `CONTRIBUTION_GUIDE.md` on how to get hold of that token (with great power...)

It's possible in the future that we can separate out the CUDA tests from the C/python tests, and automate that whole part, to automagically generate the whole "codecov" report.
