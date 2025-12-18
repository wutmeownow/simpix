# simpix

C++ starter code
* simpix_start.cpp
use make to build this example

Usage: simapix_start image1 image2 <output=out.png>

Python starter code
* simpix_start.py

Usage: simapix_start image1 image2 <output=out.png>

## instructions
Make and then run simplix with the two images of choice.<br>

Usage: ./simpix -s "srcimage" -f "tgtimage" -o "output.png" -n N -a alpha -L limit OR -t tmin <br>
Must provide paths for all three files (two in, one out), N for number of trials per pixel, alpha for annealing scheduling, limit for when to stop OR tmin for final temperature.

## results
### trump.jpg and obama.jpg -> trumpbama.pdf (1920x1080)
CMD: ./simpix -s ./trump.jpg -f ./obama.jpg -o result.png -n 1 -a 0.8 -L 0.001 <br>
Run Time: 0:03:52 <br>
Note: I tried a different command that ran for over an hour but only got slight improvement.<br>

### moon.png and milkyway.jpg -> moonyway.pdf (2000x1500)
CMD: ./simpix -s ./moon.png -f ./milkyway.jpg -o moonyway.png -n 3 -a 0.95 -L 0.00001 <br>
Run Time: 0:04:46 <br>
Note: You can see that the moon picture is not so great at replicating the photo of the milkyway, while the opposite is true. I think the milkyway pic just had much more non-black pixels while the moon pic had a dense collection of very bright white points surrounded by a sea of black, which was not great for replicating the depth of the milkyway shot.