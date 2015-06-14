set H=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361
set HB=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\bin
set HD=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\demo
set HFS=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361
set HH=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\houdini
set HHC=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\houdini\config
set HOUDINI_VERSION=14.0.361
set HSITE=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\site
set HT=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\toolkit
set HTB=C:\PROGRA~1\SIDEEF~1\HOUDIN~1.361\toolkit\bin
set Path=%HB%;%HTB%;%Path%
copy vex\include\microfacet_bsdf.h "%USERPROFILE%\Documents\houdini14.0\vex\include"
vcc --vfl-input .\vex\CVex\microfacet_sample.vfl -l .\otls\microfacet_sample.hda -A -n microfacet_sample -S microfacet_sample -N "Microfacet Sample"
vcc --vfl-input .\vex\CVex\microfacet_eval.vfl -l .\otls\microfacet_eval.hda -A -n microfacet_eval -S microfacet_eval -N "Microfacet Eval"
copy otls\microfacet_sample.hda "%USERPROFILE%\Documents\houdini14.0\otls"
copy otls\microfacet_eval.hda "%USERPROFILE%\Documents\houdini14.0\otls"