var score = new Array();
var features = new Array();
var di_score = new Array();
var di_features = new Array();
var a_fo_score = new Array();
var a_fo_features = new Array();
var g_fo_score = new Array();
var g_fo_features = new Array();
var c_fo_score = new Array();
var c_fo_features = new Array();
var u_fo_score = new Array();
var u_fo_features = new Array();


d3.csv("/static/input.csv", function(results){
	const res = results;

	const entries = Object.entries(res);

	
	//Object.entries(res).slice(33,38).map(entry => entry[1]);
for ( var i=0; i<entries.length; i++ )
{
	for ( var j=0; j<entries[i].length; j++ )
	{
	entries[i][1] = parseFloat(entries[i][1]);
	}
}
var nu = entries.slice(55,59);
for ( var i=0; i<nu.length; i++ )
{
	score.push(nu[i][1]);
	features.push(nu[i][0]);
}

var dinu = entries.slice(61,77);
for ( var i=0; i<dinu.length; i++ )
{
	di_score.push(dinu[i][1]);
	di_features.push(dinu[i][0]);
}

var basfo = entries.slice(0,32);
for ( var i=0; i<basfo.length; i++ )
{
	
	if (i <= 7){
			a_fo_score.push(basfo[i][1]);
			a_fo_features.push(basfo[i][0]);
		}
	else if (i >= 8 && i <= 15){
			g_fo_score.push(basfo[i][1]);
			g_fo_features.push(basfo[i][0]);
	}
	else if (i >= 16 && i <= 23){
			c_fo_score.push(basfo[i][1]);
			c_fo_features.push(basfo[i][0]);
	}
	else if (i >= 24 && i <= 31){
			u_fo_score.push(basfo[i][1]);
			u_fo_features.push(basfo[i][0]);
	}
}






var ctx = document.getElementById('myChart').getContext('2d');
    var chart = new Chart(ctx, {
// The type of chart we want to create
        
        type: 'bar',

        // The data for our dataset
        data: {
            labels: features,
            datasets: [{
                label: 'Nucleotide Percentage',
                backgroundColor: 'rgb(0, 99, 132)',
                borderColor: 'rgb(0, 99, 132)',
                data: score
            }]
        },


        // Configuration options go here
        options: {
        	scales: {
            yAxes: [{
                ticks: {
                    beginAtZero: true
                }
            }],
            xAxes: [{
                // Change here
            	barPercentage: 0.4
            }]
        	}
    	}
    });

 $('#0').on('click', function (e) {
        e.preventDefault;
        chart.config.data = {
            labels: di_features,
            datasets: [{
                label: "Dinucleotide Percentage",
                backgroundColor: 'rgb(19, 74, 1)',
                borderColor: 'rgb(19, 74, 1)',
                data: di_score,

            }],
        }
        chart.update();
    });
 $('#1').on('click', function (e) {
        e.preventDefault;
        chart.config.data = {
            labels: features,
            datasets: [{
                label: 'Nucleotide Percentage',
                backgroundColor: 'rgb(0, 99, 132)',
                borderColor: 'rgb(0, 99, 132)',
                data: score

            }],
        }
        chart.update();
    });
 $('#2').on('click', function (e) {
        e.preventDefault;
        chart.config.data = {
            labels: a_fo_features,
            datasets: [{
                label: 'Adenine Folding',
                backgroundColor: 'rgb(1, 11, 37)',
                borderColor: 'rgb(2, 28, 103)',
                data: a_fo_score

            }],
        }
        chart.update();
    });
 $('#3').on('click', function (e) {
        e.preventDefault;
        chart.config.data = {
            labels: g_fo_features,
            datasets: [{
                label: 'Guanine Folding',
                backgroundColor: 'rgb(1, 11, 37)',
                borderColor: 'rgb(2, 28, 103)',
                data: g_fo_score

            }],
        }
        chart.update();
    });
 $('#4').on('click', function (e) {
        e.preventDefault;
        chart.config.data = {
            labels: c_fo_features,
            datasets: [{
                label: 'Cytosine Folding',
                backgroundColor: 'rgb(1, 11, 37)',
                borderColor: 'rgb(2, 28, 103)',
                data: c_fo_score

            }],
        }
        chart.update();
    });
 $('#5').on('click', function (e) {
        e.preventDefault;
        chart.config.data = {
            labels: u_fo_features,
            
            datasets: [{
                label: 'Uracil Folding',
                backgroundColor: 'rgb(1, 11, 37)',
                borderColor: 'rgb(2, 28, 103)',
                data: u_fo_score

            }],
        }
        chart.update();
    });

});