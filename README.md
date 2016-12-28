# MSE-
Modications to the MSE Program v4.0, restructured and new input -> v4.1

These modifications were made to add "assessment frequency". In addition, some proprietary libraries were unavailable and their references have been replaced with those found in the C++ standard library. Mainly random number and statistical distribution functions found in the Rogue Wave library have been replaced.

The MSE input file now has two additional fields:

ASSESSMENT FREQUENCY

For example:

ASSESSMENT_FREQUENCY<br>
3

will analyze every three years found in the input file.

Also, in the MSR section, a field for Blim has been added after the Bmsy field:

`0.25  1  `__**`0.2`**__
