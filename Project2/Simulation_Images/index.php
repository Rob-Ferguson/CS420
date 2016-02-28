<?php

/*
	Robby Ferguson
	2/19/16
	
	This the first php script I have ever written, so please excuse
	the convolution. This script is designed to upload all of the
	simulation images produced from the CS420 Project 2 experiments
	to a basic webpage for easier comparison. Each image represents
	an individual activator/inhibitor cellular automaton produced
	from a particular set of parameters. The system was allowed to
	converge before its state space was converted to a .pgm image,
	and then I manually converted the .pgm files to .jpg images using
	an ImageMagick conversion function. For each of the three primary
	experiments run, all of the simulations are separately displayed.
	Each set of parameters was used to produce three separate AICA,
	and those are displayed together. However, if a system converged
	to 95% one state, the simulation was discarded, so those images
	are not displayed.

	To run the script, start a local server with:
		php -S localhost:8000
	Then, while that is running, navigate to http://localhost:8000
	in a browser window to see the images.
*/

print("Experiment1");
echo '<br /><br />';
$dirname = "Experiment1/";
$images = glob($dirname."*.jpg");
sort($images, SORT_NATURAL);
$prev = "500";
foreach($images as $image) {
	$parsed = explode("_", $image);
	if($prev == "500"){
		echo "Simulation $parsed[2]<br />";	
		echo '<img src="'.$image.'" border=3/>';
		print("   ");
	}
	if($parsed[2] == $prev){
		echo '<img src="'.$image.'" border=3/>';
		print("   ");
	}else{
		if($prev != "500"){
			echo "<br /><br />Simulation $parsed[2]<br />";
			echo '<img src="'.$image.'" border=3/>';
			print("   ");
		}
	}
	$prev = $parsed[2];
}

echo '<br /><br /><br />';
print("Experiment2");
echo '<br /><br />';
$dirname = "Experiment2/";
$images = glob($dirname."*.jpg");
sort($images, SORT_NATURAL);
$prev = "500";
foreach($images as $image) {
	$parsed = explode("_", $image);
	if($prev == "500"){
		echo "Simulation $parsed[2]<br />";	
		echo '<img src="'.$image.'" border=3/>';
		print("   ");
	}
	if($parsed[2] == $prev){
		echo '<img src="'.$image.'" border=3/>';
		print("   ");
	}else{
		if($prev != "500"){
			echo "<br /><br />Simulation $parsed[2]<br />";
			echo '<img src="'.$image.'" border=3/>';
			print("   ");
		}
	}
	$prev = $parsed[2];
}


echo '<br /><br /><br />';
print("Experiment3");
echo '<br /><br />';
$dirname = "Experiment3/";
$images = glob($dirname."*.jpg");
sort($images, SORT_NATURAL);
$prev = "500";
foreach($images as $image) {
	$parsed = explode("_", $image);
	if($prev == "500"){
		echo "Simulation $parsed[2]<br />";	
		echo '<img src="'.$image.'" border=3/>';
		print("   ");
	}
	if($parsed[2] == $prev){
		echo '<img src="'.$image.'" border=3/>';
		print("   ");
	}else{
		if($prev != "500"){
			echo "<br /><br />Simulation $parsed[2]<br />";
			echo '<img src="'.$image.'" border=3/>';
			print("   ");
		}
	}
	$prev = $parsed[2];
}
