<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="https://raw.githubusercontent.com/aineniamh/piranha/main/docs/piranha.svg">

    <title>${run_name} ${LANGUAGE_CONFIG["1"]}</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5zCYotlSAcp1+c8xmyTe9GYg1l9a69psu" crossorigin="anonymous">
    <script src="https://code.jquery.com/jquery-1.12.4.min.js" integrity="sha384-nvAa0+6Qg9clwYCGGPpDQLVpLNn0fRaROjHqs13t4Ggj3Ez50XnGQqc/r8MhnRDZ" crossorigin="anonymous"></script>
    <link href="https://cdn.datatables.net/1.10.25/css/jquery.dataTables.min.css" rel="stylesheet" type="text/css" />
    <link href="https://cdn.datatables.net/select/1.3.4/css/select.dataTables.min.css" rel="stylesheet" type="text/css" />
    <script src="https://cdn.datatables.net/1.10.25/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/dataTables.buttons.min.js"></script>
    <script src="https://cdn.datatables.net/select/1.3.4/js/dataTables.select.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.html5.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/buttons.print.min.js"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha384-aJ21OjlMXNL5UyIl/XNwTMqvzeRMZH2w8c5cRVpzpU8Y5bApTppSuUkhZXN0VxHd" crossorigin="anonymous"></script>
    <script src="https://cdn.jsdelivr.net/gh/rambaut/figtree.js@9880/dist/figtree.umd.js"></script>
    <script src="https://d3js.org/d3.v6.min.js"></script>
    <script src="https://sharonchoong.github.io/svg-exportJS/svg-export.min.js"></script>
    <script src="https://unpkg.com/canvg/lib/umd.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/pdfkit/js/pdfkit.min.js"></script>
    <script src="https://github.com/devongovett/blob-stream/releases/download/v0.1.3/blob-stream.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/svg-to-pdfkit/source.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/vega@5.22.1"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@5.2"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-embed@6.11.1"></script>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <% colorCodes = config["colour_map"] %>
    <% themeColor = config["colour_theme"] %>
    <style>
      body {
        padding-top: 50px;
        font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif;
      }
      table text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
      }
      header {
          display: block;
          text-align: right;
      
      }
      svg {width: 90%; height: auto;}
      .row {
          display: flex;
        }
        .column {
          padding: 10px;
          flex: 50%;
        }
      .info_box {
          background-color: #eee;
          color: #444;
          padding: 13px;
          width: 100%;
          border: none;
          text-align: left;
          outline: none;
        }
      .accordion {
          background-color: #eee;
          color: #444;
          cursor: pointer;
          padding: 13px;
          width: 100%;
          border: none;
          text-align: left;
          outline: none;
          transition: 0.4s;
        }

        .active, .accordion:hover {
          background-color: ${themeColor};
          color: white;
        }

        .accordion:after {
          content: '\002B';
          color: white;
          font-weight: bold;
          float: right;
          margin-left: 5px;
        }

        .active:after {
          content: "\2212";
        }

        .panel {
          padding: 0 13px;
          background-color: white;
          max-height: 0;
          overflow: hidden;
          transition: max-height 0.2s ease-out;
        }
      .center {
          display: block;
          margin-left: auto;
          margin-right: auto;
          width: 50%;
          }
      .svg-tooltip {
          background: rgba(69,77,93,.9);
          border-radius: .1rem;
          color: #fff;
          display: block;
          font-size: 14px;
          max-width: 320px;
          padding: .2rem .4rem;
          position: absolute;
          text-overflow: ellipsis;
          white-space: pre;
          z-index: 300;
          visibility: hidden;
    }
    .tooltip-header {
      font-size: 1.3em;
    }
    .tooltip-key {
      font-weight: bold;
    }
    .branch path{
      stroke-width:2;
      stroke: dimgrey;
      stroke-linejoin:round;
      cursor: pointer;
      }
      .branch.hovered path{
        stroke-width:4;
        stroke: dimgrey;
      }
        .node.hovered circle{
        stroke-width:5;
        stroke: dimgrey
        }
        .node text{
          font-family: "ArialNova-Light","HelveticaNeue-Light", "Helvetica Neue Light", "Helvetica Neue", Helvetica, Arial, "Lucida Grande", sans-serif; 
          font-weight: 300;
          font-size: 0.9em;
        }
      .dataTables_wrapper.no-footer .dataTables_scrollBody {
        border-top: 1px solid  rgb(148, 148, 148);
        border-bottom: none;
      }
      .svg-icon {
      display: inline-flex;
      align-self: center;
      }
      h3{
          font-size: 1em;
      }
      #toTopBtn {
      position: fixed;
      bottom: 26px;
      right: 39px;
      z-index: 98;
      padding: 21px;
      background-color: ${themeColor}
      }
      .js .cd-top--fade-out {
          opacity: .5
      }
      .js .cd-top--is-visible {
          visibility: visible;
          opacity: 1
      }
      .js .cd-top {
          visibility: hidden;
          opacity: 0;
          transition: opacity .3s, visibility .3s, background-color .3s
      }
      .cd-top {
          position: fixed;
          bottom: 20px;
          bottom: var(--cd-back-to-top-margin);
          right: 20px;
          right: var(--cd-back-to-top-margin);
          display: inline-block;
          height: 40px;
          height: var(--cd-back-to-top-size);
          width: 40px;
          width: var(--cd-back-to-top-size);
          box-shadow: 0 0 10px rgba(0, 0, 0, .05) !important;
          background: url(https://res.cloudinary.com/dxfq3iotg/image/upload/v1571057658/cd-top-arrow.svg) no-repeat center 50%;
          background-color: ${themeColor};
          background-color: hsla(var(--cd-color-3-h), var(--cd-color-3-s), var(--cd-color-3-l), 0.8)
      }
      .slidecontainer {
        width: 100%;
      }
      .colourSelect {
        background: #eee;
        border-radius: 5px;
        padding: 4px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
      }
      .slider {
        -webkit-appearance: none;
        width: 100%;
        height: 15px;
        background: #d3d3d3;
        border-radius: 5px;
        stroke: dimgrey;
        outline: none;
        opacity: 0.7;
        -webkit-transition: .2s;
        transition: opacity .2s;
      }
      .slider:hover {
        opacity: 1; 
      }
      .slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 25px;
        height: 25px;
        border-radius: 50%; 
        background: ${themeColor};
        stroke: dimgrey;
        cursor: pointer;
      }
      .slider::-moz-range-thumb {
        width: 25px;
        height: 25px;
        border-radius: 50%;
        stroke: dimgrey;
        background: ${themeColor};
        cursor: pointer;
      } 
      .tree-container{
        max-height: 1000px;
        overflow: scroll;
      }
      .label{
        display: none;
      }
      .label.show{
        display: inline;
      }
      .node.hovered .label {
          display:inline;
        }
      div.sticky {
          position: -webkit-sticky; /* Safari */
          position: sticky;
          top: 0;
        }
        .searchbar {
          border-style:solid; 
          border-color: lightgrey; 
          border-radius: 5px; 
          float:right
        }
        @media print {
                  *,
                  *::before,
                  *::after {
                    text-shadow: none !important;
                    box-shadow: none !important;
                  }
                  a:not(.btn) {
                    text-decoration: underline;
                  }
                  a[href]:after{content:none};
                  abbr[title]::after {
                    content: " (" attr(title) ")";
                  }
                  pre {
                    white-space: pre-wrap !important;
                  }
                  pre,
                  blockquote {
                    border: 1px solid #adb5bd;
                    page-break-inside: avoid;
                  }
                  tr,
                  img {
                    page-break-inside: avoid;
                  }
                  p,
                  h2,
                  h3 {
                    orphans: 3;
                    widows: 3;
                  }
                  h2,
                  h3 {
                    page-break-after: avoid;
                  }
                  .pagebreak { 
                    page-break-before: always; 
                  }
                  @page {
                    size: A4 landscape;
                    size: 287mm 210mm;
                    margin: 0.5cm;
                  }
                  body {
                    min-width: 300mm !important;
                    -webkit-print-color-adjust:exact;
                  }
                  .container {
                    min-width: 300mm !important;
                  }
                  .control {
                    display: none;
                  }
                  .searchbar {
                    display: none;
                  }
                  .badge {
                    border: 1px solid #000;
                  }
                  .table {
                    color: inherit;
                    background-color: inherit;
                    border-collapse: collapse !important;
                    display: table-row!important;
                  }
                  .table td,
                  .table th {
                    background-color: #fff !important;
                  }
                  td,td {
                      display: table-cell !important
                  }
                  .scroll-container {
                    display: none;
                }
                  .table-bordered th,
                  .table-bordered td {
                    border: 1px solid #dee2e6 !important;
                  }
                  .table-dark {
                    color: inherit;
                  }
                  .table-dark th,
                  .table-dark td,
                  .table-dark thead th,
                  .table-dark tbody + tbody {
                    border-color: #dee2e6;
                  }
                  .dataTables_scroll {
                    overflow:visible;
                  }
                  .dataTables_filter {
                    display: none;
                  }
                  .sorting_desc{
                    display: none;
                  }
                  .sorting_asc{
                    display: none;
                  }
                  .scrollX {
                    display: none;
                  }
                  .accordion {
                    display: none;
                  }
                  .panel {
                    display: none;
                  }
                }
    </style>
  </head>

  <body>
    <script>
      $(document).ready(function() {
        $(window).scroll(function() {
        if ($(this).scrollTop() > 20) {
        $('#toTopBtn').fadeIn();
        } else {
        $('#toTopBtn').fadeOut();
        }
        });
        
        $('#toTopBtn').click(function() {
        $("html, body").animate({
        scrollTop: 0
        }, 400);
        return false;
        });
        });
    </script>


    <!--Figtree.js-->
    <script type="text/javascript"> 
      const updateTableFactory = (tooltipId,metadata)=>(tipId)=>{
              const data = metadata[tipId];
              const tableDiv = d3.select(document.getElementById(tooltipId));
              //Remove table
          tableDiv.html("")
              if (data !== undefined) {
                  const visibleData = Object.keys(data).filter(d=>d!=='name');
                  tableDiv.append("h3")
                      .attr("class",'tooltip-header')
                      .text(tipId)
                      .append("hr");
                  tableDiv.selectAll("p")
                          .data(visibleData)
                          .enter()
                          .append("p")
                          .attr("class","tooltip-text")
                              .selectAll("span")
                              .data(d=>[d,data[d]])
                              .enter()
                              .append("span")
                              .attr("class",(d,i)=> i===0? "tooltip-key" :"tooltip-value")
                              .text((d,i)=>i===0? d + " : ": d);
              }
      }
      <%text>
      function addColourEventHandler(circleNodes,legend,colourSelectID,colorCodes,fig){
        d3.select(`#${colourSelectID}`).on("change", function(d){
            const selectedGroup = this.value 
            const colorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations[selectedGroup].values)
            circleNodes.attr("fill",n=>colorScale(n.annotations[selectedGroup]))

            legend.scale(colorScale)
            fig.update();
            console.log(selectedGroup);
          })
        }

      function addTraitColorEventHandler(traits,traitLegend,barSelectID,colorCodes,fig){
        
        d3.select(`#${barSelectID}`).on("change", function(d){
            const selectedGroup = this.value 
            const traitColorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations[selectedGroup].values)
            traits.attr("fill",n => traitColorScale(n.annotations[selectedGroup]))
            traitLegend.scale(traitColorScale)
            fig.update();
            console.log(selectedGroup);
          })
        }


      function addSliderEventHandler(sliderID, fig) {
          const svg = fig.svgSelection.select(function () {
              return this.parentNode;
          })
          const initialHeight = svg.attr("height");
          const maxHeight = fig.tree().externalNodes.length * 50; // 50 pixels for each tip plus a little for margins;
          if (maxHeight <= initialHeight) {
              console.log(sliderID);
              d3.select(`#${sliderID}`).select(function () {
                  return this.parentNode;
              })
                  .remove();// don't need  a slider add names
              fig.svgSelection.selectAll(".label")
                  .classed("show", true)
              return;
          }
          const heightScale = d3.scaleLinear()
              .range([initialHeight, maxHeight])
              .domain([0, 1])
          if (initialHeight / fig.tree().externalNodes.length > 12) {
              fig.svgSelection.selectAll(".label")
                  .classed("show", true)
          }
          d3.select(`#${sliderID}`).on("input", function () {
              const svgHeight = heightScale(this.value);
              //magic number!!
              svg.attr("height", svgHeight);
              fig.update();
              if (svgHeight / fig.tree().externalNodes.filter(node => !node[fig.id].ignore).length > 12) {
                  fig.svgSelection.selectAll(".label")
                      .classed("show", true)
              } else {
                  fig.svgSelection.selectAll(".label")
                      .classed("show", false)
              }
          })
      }
      </%text>
      
      function buildTree(svgID, myTreeString,tooltipID,backgroundDataString,sliderID,colourSelectID,barSelectID,colorCodes) {
          const backgroundData = JSON.parse(backgroundDataString);
          const updateTable = updateTableFactory(tooltipID, backgroundData);
          const margins = {top:50,bottom:60,left:100,right:250}
          const svg = d3.select(document.getElementById(svgID))
          svg.selectAll("g").remove();
          const nexusString = myTreeString;
          const tree = figtree.Tree.parseNexus(nexusString)[0];
          const fig = new figtree.FigTree(document.getElementById(svgID),margins, tree)
          const colorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations["source"].values)
          const traitColorScale = d3.scaleOrdinal(colorCodes).domain(fig.tree().annotations["source"].values)
          const circleNodes = figtree.circle()
                              .filter(n => !n.children)
                              .attr("r", 4)
                              .attr("fill", n => colorScale(n.annotations["source"]))
                              .hilightOnHover(10)
                              .onClick((node, i, n) => {
                                  const isSelected = d3.select(n[i]).classed("selected");
                                  fig.svgSelection.selectAll(".selected").classed("selected", false);
                                  if(isSelected){
                                      d3.select(n[i]).classed("selected", false);
                                      updateTable(null);
                                  }else{
                                      d3.select(n[i]).classed("selected", true);
                                      updateTable(node.name);
                                  }
                              });
          const legend = figtree.legend()
                                .scale(colorScale)
                                .x(-100)
                                .y(40)
          const traitLegend = figtree.legend()
                                .scale(traitColorScale)
                                .x(-150)
                                .y(40)
          const traits = figtree.traitBar()
                                .x(svg.style("width")+230)
                                .width(10)
                                .attr("fill",n => traitColorScale(n.annotations["source"]));
          fig.layout(figtree.rectangularLayout)
                  .nodes(circleNodes,
                          figtree.tipLabel(v=>v.name).attr("dx",10),
                          figtree.rectangle()
                                  .filter(n=>n[fig.id].collapsed)
                                  .attr("width",20)
                                  .attr("height",20)
                  )
                        .nodeBackgrounds(figtree.circle()
                                          .attr("r", 5)
                                .filter(n=>!n.children)
                                        )
                        .branches(figtree.branch()
                                    .hilightOnHover(20) 
                                    .collapseOnClick()
                                    .on("click",()=>{
                                      const svgHeight = fig.svgSelection.select(function() { return this.parentNode; }).attr("height");
                                      if(svgHeight/fig.tree().externalNodes.filter(node=>!node[fig.id].ignore).length>12){
                                        fig.svgSelection.selectAll(".label")
                                          .classed("show",true)
                                      }else{
                                        fig.svgSelection.selectAll(".label")
                                        .classed("show",false)
                                      }
                                    })
                            )
                            .feature(
                                    figtree.scaleBar()
                                      .direction("x")
                                      .length(1/900)
                                      .x(-60)
                                      .y(-30)
                                      // .y(fig.svgSelection.select(function() { return this.parentNode; }).attr("height")-margins.top-margins.bottom+20)
                                      .title({text:"~1 SNP",
                                      yPadding:10})
                                        )
                            .feature(legend)
                            // .feature(traits)
        addSliderEventHandler(sliderID,fig);
        addColourEventHandler(circleNodes,legend,colourSelectID,colorCodes,fig);
        // addTraitColorEventHandler(traits,traitLegend,barSelectID,colorCodes,fig)
      }
    </script>



<script>
  function exportImageSVG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadSvg(document.querySelector(svgID), name);
      };
  };
  function exportImagePNG(buttonID,svgID,name){
      document.querySelector(buttonID).onclick = function(){
          svgExport.downloadPng(document.querySelector(svgID), name);
      };
  };
</script>
    <div class="container">
      <a href="#" id="toTopBtn" class="cd-top text-replace js-cd-top cd-top--is-visible cd-top--fade-out" data-abc="true"></a>
      <div>
      <header class="piranha-header">
        <div class="col-sm-4" style="text-align: left;"><a href="https://polionanopore.org/">
          <img class="piranha-logo" src="https://raw.githubusercontent.com/aineniamh/piranha/main/docs/poseco.svg" vertical-align="left" width="30" height="30"></img></a>
          PoSeCo | <small class="text-muted">${LANGUAGE_CONFIG["2"]}</small>
        </div>
        <div class="col-sm-8" style="text-align: right;"><a style="color:#000" href="https://github.com/polio-nanopore/piranha/">
          <strong style="color:#000"> piranha </strong> | <small class="text-muted">${LANGUAGE_CONFIG["3"]}</small></a>
        </div>
        <br>
        <hr>
      </header>
        
      <h1>${LANGUAGE_CONFIG["4"]} ${run_name} <small class="text-muted" style="color:${themeColor}">${date}</small></h1>
      <br>
      %if config["username"]!="":
        <h3><strong>${LANGUAGE_CONFIG["5"]}</strong> | ${config["username"].lstrip("'").rstrip("'")}</h3>
      %endif
      %if config["institute"]!="":
        <h3><strong>${LANGUAGE_CONFIG["6"]}</strong> | ${config["institute"].lstrip("'").rstrip("'")}</h3>
      %endif
      %if config["notes"]!="":
        <hr>
        <p>${config["notes"].lstrip("'").rstrip("'")}</p>
        <hr>
      %endif
      <div class="info_box">
        <p>${LANGUAGE_CONFIG["52"]}</p>
      </div>
      <br>
      %if flagged_high_npev and config["sample_type"]=="environmental":
      <h3><strong>${LANGUAGE_CONFIG["7"]}</strong> | ${LANGUAGE_CONFIG["8"]}</h3>
      <p style="white-space:wrap; word-wrap:break-word; overflow:scroll; border-width:2px; border-style:solid; border-color:#e68781; padding: 1em;">
      %for sample in flagged_high_npev:
        - ${sample}<br>
      %endfor
    </p>
    %endif  
      <br>
      <h3><strong>${LANGUAGE_CONFIG["9"]} 1</strong> | ${LANGUAGE_CONFIG["10"]} </h3>
      <button class="accordion">${LANGUAGE_CONFIG["11"]}</button>
        <div class="panel">
          <div class="row">
            <div class="col-sm-2" ><strong>${LANGUAGE_CONFIG["11"]}: </strong></div>
            <div class="col-sm-4" id="tableExportID1"></div>
            <div class="col-sm-3" ><strong>${LANGUAGE_CONFIG["13"]}: </strong></div>
            <div class="col-sm-2" ><button><a href="${detailed_csv_out}" download="detailed_run_report.csv">${LANGUAGE_CONFIG["12"]}</a></button></div>
          </div>
        </div>
        <table class="display nowrap" id="myTable1">
          <thead>
            <tr>
             <% header = LANGUAGE_CONFIG["14"] %>
              %for col in header:
                <th>${col.title().replace("_"," ")}</th>
              %endfor
              <th>${LANGUAGE_CONFIG["15"]} (${config["analysis_mode"].upper()})</th>
            </tr>
          </thead>
          <tbody>
            <% import collections %>
            <% summary_bcodes = collections.defaultdict(list) %>
            % for row in data_for_report["summary_table"]:
                %for col in config["summary_table_header"]:

                  %if col=="sample":
                    <% this_barcode = row["barcode"] %>
                    <% this_reference = row["reference_group"] %>
                    <% summary_bcodes[this_barcode].append(this_reference) %>
                    <td><a href="./barcode_reports/${this_barcode}_report.html" target="_blank" style="color:${themeColor}"><strong>${row[col]}</strong></a></td>
                  %else:
                    <td>${row[col]}</td>
                  %endif
                %endfor
                <td><a download href="published_data/${this_barcode}/${this_barcode}.consensus.fasta" style="color:${themeColor}"><strong>${LANGUAGE_CONFIG["16"]}</strong></a></td>
              </tr>
            % endfor
          </tbody>
        </table>
        <script type="text/javascript">
          $(document).ready( function () {
              var table = $('#myTable1').DataTable({
                select: {
                        style: 'multi'
                    },
                'iDisplayLength': 100,
                "paging": false,
                "border-bottom":false,
                "bInfo" : false,
                dom: 'frtip',
                buttons: ["copy","csv","print"]
              });
              table.buttons().container().appendTo( $('#tableExportID1') );
              
            } );
        </script>
      <div class="pagebreak"> </div>
      <h3><strong>${LANGUAGE_CONFIG["9"]} 2</strong> | ${LANGUAGE_CONFIG["17"]} </h3>
      <button class="accordion">${LANGUAGE_CONFIG["11"]}</button>
        <div class="panel">
          <div class="row">
            <div class="col-sm-2" ><strong>${LANGUAGE_CONFIG["11"]}: </strong></div>
            <div class="col-sm-8" id="tableExportID2"></div>
          </div>
        </div>
        <table class="display nowrap" id="myTable2">
          <thead>
            <tr>
              %for col in config["composition_table_header"]:
                %if col=="sample":
                <th style="width:10%;">${LANGUAGE_CONFIG["18"]}</th>
                %elif col=="barcode":
                <th style="width:10%;">${LANGUAGE_CONFIG["58"]}</th>
                %elif col=="unmapped":
                <th style="width:10%;">${LANGUAGE_CONFIG["19"]}</th>
                %else:
                <th style="width:10%;">${col.replace("_"," ")}</th>
                %endif
              %endfor
            </tr>
          </thead>
          <tbody>
            
            % for row in data_for_report["composition_table"]:
              %if row["sample"] not in config["negative_control"] and row["sample"] not in config["positive_control"]:
              <% this_barcode = row["barcode"] %>
              
              <tr>
                %for col in config["composition_table_header"]:

                  %if col=="sample":
                    
                    %if this_barcode in summary_bcodes:
                    <td><a href="./barcode_reports/${this_barcode}_report.html" target="_blank" style="color:${themeColor}">${row[col]}</a></td>
                    %else:
                    <td style="color:${themeColor}">${row[col]}</td>
                    %endif

                  %elif col!="barcode":

                    %if this_barcode in summary_bcodes:
                      <% these_references = summary_bcodes[this_barcode] %>
                        %if col in these_references:
                          <td style="color:${themeColor}"><strong>${row[col]}</strong></td>
                        %else:
                          <td>${row[col]}</td>
                        %endif
                    %else:
                      <td>${row[col]}</td>
                    %endif

                  %else:
                    <td>${row[col]}</td>
                  %endif

                %endfor
              </tr>
              %endif
            % endfor
          </tbody>
        </table>
        <br>
        <p>*${LANGUAGE_CONFIG["57"]} ${config["composition_not_detected"]}</p>
        <script type="text/javascript">
          $(document).ready( function () {
              var table = $('#myTable2').DataTable({
                'iDisplayLength': 100,
                "paging": false,
                "border-bottom":false,
                "bInfo" : false,
                dom: 'frtip',
                buttons: ["copy","csv","print"]
              });
              table.buttons().container().appendTo( $('#tableExportID2') );
              
            } );
        </script>

    <div class="pagebreak"> </div>
    %if flagged_seqs:
    <h3><strong>${LANGUAGE_CONFIG["9"]} 3</strong> | ${LANGUAGE_CONFIG["20"]} </h3>
    <button class="accordion">${LANGUAGE_CONFIG["11"]}</button>
      <div class="panel">
        <div class="row">
          <div class="col-sm-2" ><strong>${LANGUAGE_CONFIG["11"]}: </strong></div>
          <div class="col-sm-8" id="tableExportID3"></div>
        </div>
      </div>
      
      <table class="display nowrap" id="myTable3">
        <thead>
          <tr>
              <th style="width:30%;">${LANGUAGE_CONFIG["21"]}</th>
              <th style="width:70%;">${LANGUAGE_CONFIG["22"]}</th>
          </tr>
        </thead>
        <tbody>
          <% set_count = 0 %>
          % for set in flagged_seqs:
            <% set_count += 1 %>
            <tr>
              <td>${set_count}</td>
              <td style="overflow:scroll;">
                % for seqid in set:
                  ${seqid}<br>
                %endfor
              </td>
            </tr>

          % endfor
        </tbody>
      </table>
      <br>
      <script type="text/javascript">
        $(document).ready( function () {
            var table = $('#myTable3').DataTable({
              'iDisplayLength': 100,
              "paging": false,
              "border-bottom":false,
              "bInfo" : false,
              dom: 'frtip',
              buttons: ["copy","csv","print"]
            });
            table.buttons().container().appendTo( $('#tableExportID3') );
            
          } );
      </script>
      %endif
      % if show_control_table:
        <div class="pagebreak"> </div>
        <h3><strong>${LANGUAGE_CONFIG["9"]} 4</strong> | ${LANGUAGE_CONFIG["23"]} </h3>
          <table class="table">
            <thead class="thead-light">
              <tr>
                <th style="width:8%;">${LANGUAGE_CONFIG["24"]}</th>
                %for col in config["composition_table_header"]:
                    %if col=="sample":
                    <th>${LANGUAGE_CONFIG["18"]}</th>
                    %elif col=="barcode":
                    <th>${LANGUAGE_CONFIG["58"]}</th>
                    %elif col=="unmapped":
                    <th>${LANGUAGE_CONFIG["19"]}</th>
                    %else:
                    <th>${col.replace("_"," ")}</th>
                    %endif
                %endfor
              </tr>
            </thead>
            <tbody>
              % for row in data_for_report["composition_table"]:
                %if row["sample"] in data_for_report["control_status"]:
                  <% control_status = data_for_report["control_status"][row["sample"]] %>
                  %if control_status:
                    <tr style="background-color:rgba(25, 67, 76, 0.3)">
                    <td>
                      <span class="glyphicon">&#xe013;</span>
                    </td>
                  %else:
                    <tr style="background-color:rgba(230, 135, 129, 0.3)">
                      <td></td>
                  %endif
                      %for col in config["composition_table_header"]:
                        %if col not in ["barcode","sample"]:
                          %if int(row[col])>config["min_read_depth"]:
                            <td><strong>${row[col]}</strong></td>
                          %else:
                            <td>${row[col]}</td>
                          %endif
                        %else:
                          <td>${row[col]}</td>
                        %endif
                        
                      %endfor
                    </tr>
                %endif
              %endfor
            </tbody>
          </table>
      %endif
      </div>
    <br>
    
  <div class="pagebreak"> </div>
  <div id="plateViz"></div>
    <script>
       var vlSpec_plate = {
         "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
         "width": 450,
         "height": 350,
         "data": {"values": ${plate_json}},
       "params": [
           {"name":"filterBy",
             "value":"All",
             "bind":{"input":"select",
                     "options":${positive_types},
                     "labels":${positive_types},
                     "name":"Filter by: "
                    }}
         ],
       "layer":[
             {
               "transform": [{"calculate": "datum[filterBy]", "as": "Reads present"}],
               "mark": {"type":"circle","size":500},
               "encoding": {
                 "x": {"field": "x", "type": "nominal",
                 "title": "",
                         "axis": {"grid": false,
                                 "labelFont":"Helvetica Neue",
                                 "labelFontSize":18,
                                 "labelAngle": 0
                                 }},
                 "y": {"field": "y", "type": "ordinal",
                 "title": "",
                         "axis": {"grid": false,
                                 "labelFont":"Helvetica Neue",
                                 "labelFontSize":18
                                 }},
                 "fill": {"field": "EV reads present",
                          "scale": {"range": ["#e68781", "#48818d", "#b2b2b2"]},
                          "sort": ["Present","Absent","N/A"]
                        },
               "tooltip": [
                 {"field": "Barcode", "type": "nominal"},
                 % for t in positive_types:
                   {"field": "${t}", "type": "nominal"},
                 % endfor
               ]
               }
             }]
         };
               vegaEmbed('#plateViz', vlSpec_plate, {renderer: "svg"})
                     .then(result => console.log(result))
                     .catch(console.warn);
     </script>
     
     <h3><strong>${LANGUAGE_CONFIG["25"]} 1</strong> | ${LANGUAGE_CONFIG["26"]}</h3>
     <br>
    
      
      %if config["run_phylo"]:
      <% figure_count = 1 %>
      %for reference_group in phylo_data:
      <% print(reference_group) %>
      <% figure_count +=1 %>
      <div class="pagebreak"> </div>

      <button class="accordion">${LANGUAGE_CONFIG["27"]}</button>
        <div class="panel">
          <div class="row">
            <div class="column">
              <div class="slider-block" id="slider_${reference_group}">
                <p>${LANGUAGE_CONFIG["28"]}</p>
                <input class="slider" type="range" id="rangeinput_${reference_group}"  min="0" max="1" style="width: 100px" step="0.01" value="0" />
                <span class="highlight"></span>
              </div>
            </div>
            <div class="column">
              <div>
              <p>${LANGUAGE_CONFIG["29"]}</p>
              <select class="colourSelect" id="colourSelect_${reference_group}">
                % for annotation in config["tree_annotations"].rstrip(" ").split(" "):
                  <option value="${annotation}">${annotation.title()}</option>
                % endfor
              </select>
              </div>
            </div>
          </div>
        </div>
      <div class="row tree-container">
        <div class="col-xs-7">
          <svg class="tree_svg" width="700" height="400" id="tree_${reference_group}"></svg>
        </div>
        <div class="col-xs-4 sticky" id="tooltip_${reference_group}">
        </div> 
        
        <script type="text/javascript">
          buildTree("tree_${reference_group}", 
                    `${phylo_data[reference_group]['nexus']}`,
                    `tooltip_${reference_group}`,
                    `${background_data}`,
                    "rangeinput_${reference_group}",
                    "colourSelect_${reference_group}",
                    "barSelect_${reference_group}",
                    ${colorCodes});
        </script> 
      </div> 
      <h3><strong>${LANGUAGE_CONFIG["25"]} ${figure_count}</strong> | ${reference_group} ${LANGUAGE_CONFIG["30"]}</h3>
      <hr>
    %endfor
    %endif

    <div class="pagebreak"> </div>
    <br>
    %if show_control_table:
      <h3><strong>${LANGUAGE_CONFIG["9"]} 5</strong> | ${LANGUAGE_CONFIG["31"]} </h3>
    %else:
    <h3><strong>${LANGUAGE_CONFIG["9"]} 4</strong> | ${LANGUAGE_CONFIG["31"]} </h3>
    %endif
    <button class="accordion">${LANGUAGE_CONFIG["11"]}</button>
      <div class="panel">
        <div class="row">
          <div class="col-sm-2" ><strong>${LANGUAGE_CONFIG["11"]}: </strong></div>
          <div class="col-sm-8" id="tableExportID5"></div>
        </div>
      </div>
        <table class="display nowrap" id="myTable5">
          <thead>
            <tr>
              <th style="width:30%;">${LANGUAGE_CONFIG["32"]}</th>
              <th style="width:70%;">${LANGUAGE_CONFIG["33"]}</th>
            </tr>
          </thead>
          <tbody>
            % for item in config["configuration_table_fields"]:
              <tr>
                <td><strong>${item}</strong></td>
                %if type(config[item]) == list:
                <% new_item = ", ".join(config[item]) %>
                <td>${new_item}</td>
                %else:
                <td>${config[item]}</td>
                %endif
            %endfor
          </tbody>
        </table>
        <script type="text/javascript">
          $(document).ready( function () {
              var table = $('#myTable5').DataTable({
                select: {
                        style: 'multi'
                    },
                'iDisplayLength': 100,
                "paging": false,
                "border-bottom":false,
                "bInfo" : false,
                dom: 'frtip',
                buttons: ["copy","csv","print"]
              });
              table.buttons().container().appendTo( $('#tableExportID5') );
              
            } );
        </script>
      <br>
      <script>
        var acc = document.getElementsByClassName("accordion");
        var i;
        for (i = 0; i < acc.length; i++) {
              acc[i].addEventListener("click", function() {
                this.classList.toggle("active");
                var panel = this.nextElementSibling;
                if (panel.style.maxHeight) {
                  panel.style.maxHeight = null;
                } else {
                  panel.style.maxHeight = panel.scrollHeight*1.2 + "px";
                } 
              });
            }
      </script>
    <div id="citation_box" class="info_box">
      <p>${LANGUAGE_CONFIG["53"]}</p>
      <p>
        <strong>O’Toole Á, Colquhoun R, Ansley C, Troman C, Maloney D, Vance Z, Akello J, Bujaki E, Majumdar M, Khurshid A, Arshad Y, Alam MM, Martin J, Shaw A, Grassly N, Rambaut A</strong> (2023) Automated detection and classification of polioviruses from nanopore sequencing reads using piranha. <i>Virus Evolution</i> <a style='color:#e68781' href="https://doi.org/10.1093/ve/veae023">https://doi.org/10.1093/ve/veae023</a>
      </p>
    </div>
    <footer class="page-footer">
      <div class="container-fluid text-right text-md-right">
        <hr>
        <div class="row">
          <div class="col-sm-1">
            <p>
            <img class="piranha-logo" src="https://raw.githubusercontent.com/aineniamh/piranha/main/docs/piranha.svg?token=GHSAT0AAAAAABHOJJPRFXUJULWKXQA5SVOMYP46NBA" vertical-align="left" width="50" height="50"></img>
            <p>
          </div>
          <div class="col-sm-11" style="text-align: right;">
            piranha ${version} | <small class="text-muted">${LANGUAGE_CONFIG["3"]}</small> <br><small class="text-muted">GNU General Public License v3.0</small>
          </div>
        <br><br>
        </p>
      </div>
    </footer>
  </div>
</body>
</html>
