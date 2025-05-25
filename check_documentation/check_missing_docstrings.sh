pylint --load-plugins=pylint.extensions.docparams\
       --disable=all \
       --enable=missing-docstring,missing-param-doc,missing-type-doc,missing-return-doc \
        ../wodenpy ../scripts/run_woden.py > docstring_issues.txt