var title = "Dario Webserver"; 
function addToFavorites(url) { 
	if (window.sidebar) { // Mozilla Firefox Bookmark
		window.sidebar.addPanel(title, url,"");
	} else if( window.external ) { // IE Favorite
		window.external.AddFavorite( url, title); }
	else if(window.opera && window.print) { // Opera Hotlist
		return true; }
  // if (window.external) { window.external.AddFavorite(url,title)} 
  // else { alert("Sorry! Your browser doesn't support this feature.");} 
} 

// Function to return if an array contains an element or not
function contains(a, obj) {
  for(var i = 0; i < a.length; i++) {
    if(a[i] === obj){
      return true;
    }
  }
  return false;
}

function button01() {
  // Print error if test data is not available
  if (typeof(species_with_testdata) != "undefined" && document.mtdb.use_test_data == true) {
    // alert(document.mtdb.code.value);
    if (!contains(species_with_testdata, document.mtdb.code.value)) {
      alert("There is currently no test data available for the selected species.");
      return 0;
    };
  };

  if(document.mtdb.coverage_file.value != "" || document.mtdb.use_test_data.checked == true) {
    document.mtdb.submit();
  }
  else {
    alert("You must submit a file!");
  }
}
