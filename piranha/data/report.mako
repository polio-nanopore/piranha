<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">

    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="https://raw.githubusercontent.com/aineniamh/piranha/main/docs/piranha.svg">

    <title>${run_name} report</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5zCYotlSAcp1+c8xmyTe9GYg1l9a69psu" crossorigin="anonymous">
    <script src="https://code.jquery.com/jquery-1.12.4.min.js" integrity="sha384-nvAa0+6Qg9clwYCGGPpDQLVpLNn0fRaROjHqs13t4Ggj3Ez50XnGQqc/r8MhnRDZ" crossorigin="anonymous"></script>
    <link href="https://cdn.datatables.net/1.10.25/css/jquery.dataTables.min.css" rel="stylesheet" type="text/css" />
    <script src="https://cdn.datatables.net/1.10.25/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.7.1/js/dataTables.buttons.min.js"></script>
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
    <script src="https://cdn.jsdelivr.net/npm/vega@5.16.0"></script>
  <script src="https://cdn.jsdelivr.net/npm/vega-lite@4.15.0"></script>
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
        <div class="col-sm-4" style="text-align: left;">
          <img class="piranha-logo" src="https://raw.githubusercontent.com/aineniamh/piranha/main/docs/poseco.svg" vertical-align="left" width="30" height="30"></img>
          PoSeCo | <small class="text-muted">Poliovirus Sequencing Consortium</small>
        </div>
        <div class="col-sm-8" style="text-align: right;">
          piranha | <small class="text-muted">Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis</small>
        </div>
        <br>
        <hr>
      </header>
        
      <h1>${run_name} report <small class="text-muted" style="color:${themeColor}">${date}</small></h1>
      %if config["username"]!="":
        <h3><strong>User</strong> | ${config["username"].lstrip("'").rstrip("'")}</h3>
      %endif
      %if config["institute"]!="":
        <h3><strong>Institute</strong> | ${config["institute"].lstrip("'").rstrip("'")}</h3>
      %endif
      <br>
      <h3><strong>Table 1</strong> | Sample summary information </h3>
      <button class="accordion">Export table</button>
        <div class="panel">
          <div class="row">
            <div class="col-sm-2" ><strong>Export table: </strong></div>
            <div class="col-sm-8" id="tableExportID1"></div>
          </div>
        </div>
        <table class="display nowrap" id="myTable1">
          <thead>
            <tr>
              %for col in config["summary_table_header"]:
                <th>${col.title().replace("_"," ")}</th>
              %endfor
            </tr>
          </thead>
          <tbody>
            % for row in data_for_report["summary_table"]:
                %for col in config["summary_table_header"]:

                  %if col=="sample":
                    <% this_barcode = row["barcode"] %>
                    <td><a href="./barcode_reports/${this_barcode}_report.html" target="_blank" style="color:${themeColor}">${row[col]}</a></td>
                  %else:
                    <td>${row[col]}</td>
                  %endif

                %endfor
              </tr>
            % endfor
          </tbody>
        </table>
        <script type="text/javascript">
          $(document).ready( function () {
              var table = $('#myTable1').DataTable({
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
      <h3><strong>Table 2</strong> | Composition of samples </h3>
      <button class="accordion">Export table</button>
        <div class="panel">
          <div class="row">
            <div class="col-sm-2" ><strong>Export table: </strong></div>
            <div class="col-sm-8" id="tableExportID2"></div>
          </div>
        </div>
        <table class="display nowrap" id="myTable2">
          <thead>
            <tr>
              %for col in config["composition_table_header"]:
                <th style="width:10%;">${col.title().replace("_"," ")}</th>
              %endfor
            </tr>
          </thead>
          <tbody>
            % for row in data_for_report["composition_table"]:
              %if row["sample"] not in [config["negative"],config["positive"]]:
              <tr>
                %for col in config["composition_table_header"]:

                  %if col=="sample":
                    <% this_barcode = row["barcode"] %>
                    <td><a href="./barcode_reports/${this_barcode}_report.html" target="_blank" style="color:${themeColor}">${row[col]}</a></td>
                  %elif col!="barcode":

                    %if int(row[col])>config["min_read_depth"]:
                      <td style="color:${themeColor}">${row[col]}</td>
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
      % if show_control_table:
        <div class="pagebreak"> </div>
        <h3><strong>Table 3</strong> | Controls </h3>
          <table class="table">
            <thead class="thead-light">
              <tr>
                <th>Pass</th>
                %for col in config["composition_table_header"]:
                  <th>${col.title().replace("_"," ")}</th>
                %endfor
              </tr>
            </thead>
            <tbody>
              % for row in data_for_report["composition_table"]:
                %if row["sample"] in [config["negative"],config["positive"]]:
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
            piranha ${version} | <small class="text-muted">Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis</small> <br><small class="text-muted">GNU General Public License v3.0</small>
          </div>
        <br><br>
        </p>
      </div>
    </footer>
  </div>
</body>
</html>
