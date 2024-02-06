
import { select } from 'd3';
import { showKOTable } from './init';

export function addColumnPicker(row, columns, selectedColumn) {
  const grp = row.append('div').attr('class', 'col-lg-2 form-group columnPicker');
  grp.append('label').text('Column');
  grp.append('select')
    .attr('class', 'form-control')
    .on('change', function changeColumn() {
      const body = select('#main');
      const data = d[this.selectedIndex];
      document.getElementById("ItemPreview").src = "data:image/png;base64," + data.image;
      showKOTable(body, data);
    })
    .selectAll('option')
    .data(columns)
    .enter()
      .append('option')
      .attr('value', d => d)
      .text(d => d)
      .property('selected', d => (d === selectedColumn));
  return grp;
}

export function addDownloadLinks(sel, svg) {
}
