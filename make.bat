set HOUDINI_VERSION=14.0.361
set HROOT=C:\Program Files\Side Effects Software
set H=%HROOT%\Houdini %HOUDINI_VERSION%
set HB=%H%\bin
set HD=%H%\demo
set HFS=%H%
set HH=%HFS%\houdini
set HHC=%HH%\config
set HSITE=%HFS%\site
set HT=%HFS%\toolkit
set HTB=%HT%\bin
set Path=%HB%;%HTB%;%Path%
rem copy /Y vex\include\microfacet_bsdf.h "%USERPROFILE%\Documents\houdini14.0\vex\include"
vcc -I vex\include --vfl-input .\vex\CVex\microfacet_sample.vfl -l .\otls\microfacet_sample.hda -A -n microfacet_sample -S microfacet_sample -N "Microfacet Sample"
vcc -I vex\include --vfl-input .\vex\CVex\microfacet_eval.vfl -l .\otls\microfacet_eval.hda -A -n microfacet_eval -S microfacet_eval -N "Microfacet Eval"
