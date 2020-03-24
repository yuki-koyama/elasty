clang-format src/*.cpp include/elasty/*.hpp tests/*.*pp examples/**/*.*pp --output-replacements-xml | grep -c "<replacement " >/dev/null
if [ $? -ne 1 ]; then echo "Ill-styled code detected" && return 1; else return 0; fi
