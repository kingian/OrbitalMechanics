
//self calling anonymous function
(function ()
{
	var a = {}
	
	function addto(a,b)
	{
		return a + b;	
	}
	
	
	x = addto(1,2)
}())