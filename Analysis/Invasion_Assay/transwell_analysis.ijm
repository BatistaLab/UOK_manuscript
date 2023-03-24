dir1 = getDirectory("/Users/fitzsimmonscm/Documents/Phenotype Assays/Invasion Assays/Transwell Invasion/invaded_cells/Cropped");
dir2 = getDirectory("/Users/fitzsimmonscm/Documents/Phenotype Assays/Invasion Assays/Transwell Invasion/invaded_cells/Cropped/Results");
list = getFileList(dir1);
setBatchMode(true);
for (i=0; i<list.length; i++) {

showProgress(i+1, list.length);
open(dir1+list[i]);
title = getTitle();
//run("Threshold...");
setThreshold(0, 26981);
setOption("BlackBackground", true);
run("Convert to Mask", "method=Default background=Dark black");
run("Watershed", "slice");
run("Analyze Particles...", "size=30-100 show=[Bare Outlines] display clear summarize slice");
selectWindow("Results");
saveAs("Results", dir2 + title +".csv");
saveAs("Tiff", dir2 + title +".tif");
 close();

}
selectWindow("Summary");
saveAs("Text", dir2 + "summary.txt");



