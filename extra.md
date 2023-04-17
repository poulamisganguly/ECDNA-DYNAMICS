if you cannot checkout branch then force it :
`git checkout -f biasedDoubling'

if you get weird old mode new mode files in `git status' then run:
`git config core.filemode false'
Check [this](https://stackoverflow.com/questions/1257592/how-do-i-remove-files-saying-old-mode-100755-new-mode-100644-from-unstaged-cha)

on my work computer to make the executable I have to run:
make CXX=/usr/lib64/ccache/x86_64-redhat-linux-g++
