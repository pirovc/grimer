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
    //If numeric, replace NaN with sortable value (-Infinity) to be at the end of the sorted array
	if (sort_mode=="numeric"){
		sort_col = sort_col.map(function(v){ return isNaN(v) ? -Infinity : v })
	}

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

function table_to_tsv(source, cols, headers, selected) {

	var rows_idx = []
	if(selected==true){
		//remove undefined from selected if present
		rows_idx = source.selected.indices.filter(function( element ) {
		   return element !== undefined;
		});
	}
	else{
		// include all rows
		for (let i = 0; i < source.get_length(); i++) {
			rows_idx.push(i);
		}
	}

	const lines = [headers.join('\t')]
    for (let i = 0; i < rows_idx.length; i++) {
    	let row = [];
        for (let j = 0; j < cols.length; j++) {
            row.push(source.data[cols[j]][rows_idx[i]].toString())
        }
        lines.push(row.join('\t'))
    }
    return lines.join('\n').concat('\n')
}

function save_file(filename, filetext){
    const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })
    //addresses IE
    if (navigator.msSaveBlob) {
        navigator.msSaveBlob(blob, filename)
    } else {
        const link = document.createElement('a')
        link.href = URL.createObjectURL(blob)
        link.download = filename
        link.target = '_blank'
        link.style.visibility = 'hidden'
        link.dispatchEvent(new MouseEvent('click'))
    }
}