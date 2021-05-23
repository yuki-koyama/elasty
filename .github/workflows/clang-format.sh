clang-format src/*.cpp include/elasty/*.hpp tests/*.*pp examples/**/*.*pp --output-replacements-xml
clang-format src/*.cpp include/elasty/*.hpp tests/*.*pp examples/**/*.*pp --output-replacements-xml | grep -c "<replacement " >/dev/null
if [ $? -ne 1 ]; then echo "Ill-styled code detected" && exit 1; fi
