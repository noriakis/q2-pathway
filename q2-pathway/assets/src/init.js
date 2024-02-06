
/* global d */
import { select } from 'd3';
import { addDownloadLinks, addColumnPicker } from './toolbar';


export default function init(groupIndex) {
  const data = d[groupIndex];
  const columns = d.map(d => d.column);

  // DOM
  const body = select('#main');
  const plotRow = body.insert('div', ':first-child').attr('class', 'viz row');
  const plotDiv = plotRow.append('div').attr('class', 'col-lg-12');
  const controlsRow = plotDiv.append('div').attr('class', 'controls row');
  const img = plotRow.append("center").append("img").attr('id', "ItemPreview").attr("src", "");

  body.insert('h1', ':first-child').text('KEGG PATHWAY image');

  // Plot main
  document.getElementById("ItemPreview").src = "data:image/png;base64," + data.image;

  // CONTROLS
//  addDownloadLinks(controlsRow, svg);
  addColumnPicker(controlsRow, columns, data.column);

  // Show KO Table
  showKOTable(body, data);
}

export function showKOTable(body, data) {
  select('#ko-table').html(data.table);
  select('#ko-csv').attr('href', data.KOCSVPath);
}
