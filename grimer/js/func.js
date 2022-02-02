function sort_numeric(a, b){ return a - b; }
function sort_string(a, b){ return a.localeCompare(b); }

function grimer_sort(factors, sort_col, sort_mode="numeric", desc=false, group_col1=[], group_col2=[], index=[]) {
    //sort_mode : numeric, string

    // subset data if index provided
    if(index.length){
    	factors = index.map( s => factors[s] );
		sort_col = index.map( s => sort_col[s] );
		if(group_col1.length){
			group_col1 = index.map( s => group_col1[s] );
		}
		if(group_col2.length){
			group_col2 = index.map( s => group_col2[s] );
		}
    }

    // Generate numerical index to sort arrays
    var idx = new Array(factors.length);
    for (var i = 0; i < idx.length; ++i) idx[i] = i;
    //If numeric, replace NaN with sortable value (false)
	if (sort_mode=="numeric") sort_col = sort_col.map(function(v){ return isNaN(v) ? false : v })

	if(group_col1.length && group_col2.length){
		if (sort_mode=="numeric" && desc==false)
			idx.sort((a, b) => sort_string(group_col2[a],group_col2[b]) || sort_string(group_col1[a],group_col1[b]) || sort_numeric(sort_col[b],sort_col[a]));
		else if (sort_mode=="numeric" && desc==true)
			idx.sort((a, b) => sort_string(group_col2[a],group_col2[b]) || sort_string(group_col1[a],group_col1[b]) || sort_numeric(sort_col[a],sort_col[b]));
		else if (sort_mode=="string" && desc==false)
			idx.sort((a, b) => sort_string(group_col2[a],group_col2[b]) || sort_string(group_col1[a],group_col1[b]) || sort_string(sort_col[a],sort_col[b]));
		else if (sort_mode=="string" && desc==true)
			idx.sort((a, b) => sort_string(group_col2[a],group_col2[b]) || sort_string(group_col1[a],group_col1[b]) || sort_string(sort_col[b],sort_col[a]));
	}else if(group_col1.length){
		if (sort_mode=="numeric" && desc==false)
			idx.sort((a, b) => sort_string(group_col1[a],group_col1[b]) || sort_numeric(sort_col[b],sort_col[a]));
		else if (sort_mode=="numeric" && desc==true)
			idx.sort((a, b) => sort_string(group_col1[a],group_col1[b]) || sort_numeric(sort_col[a],sort_col[b]));
		else if (sort_mode=="string" && desc==false)
			idx.sort((a, b) => sort_string(group_col1[a],group_col1[b]) || sort_string(sort_col[a],sort_col[b]));
		else if (sort_mode=="string" && desc==true)
			idx.sort((a, b) => sort_string(group_col1[a],group_col1[b]) || sort_string(sort_col[b],sort_col[a]));
	}else{
		if (sort_mode=="numeric" && desc==false)
			idx.sort((a, b) => sort_numeric(sort_col[b],sort_col[a]));
		else if (sort_mode=="numeric" && desc==true)
			idx.sort((a, b) => sort_numeric(sort_col[a],sort_col[b]));
		else if (sort_mode=="string" && desc==false)
			idx.sort((a, b) => sort_string(sort_col[a],sort_col[b]));
		else if (sort_mode=="string" && desc==true)
			idx.sort((a, b) => sort_string(sort_col[b],sort_col[a]));
	}
	
    var sorted_factors = new Array(idx.length);
    for (var i = 0; i < idx.length; ++i) sorted_factors[i] = factors[idx[i]];
    return sorted_factors;
}
