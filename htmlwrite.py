# a document with assorted functions for creating the exoplanet webpage
import urllib
import os

# function to convert degrees into ra hh:mm:ss string
def _deg_to_ra_hms(angle_deg, delim = ':', round_seconds = False):
    hours = (angle_deg - (angle_deg % 15)) / 15
    minutes = (angle_deg / 15 - hours) * 60
    seconds = (minutes % 1) * 60
     
    # decide to round the seconds (for display) or leave unrounded (internal)
    if round_seconds: return str(int(hours)) + delim + str(int(minutes)) + delim + str(round(seconds))
    else: return str(int(hours)) + delim + str(int(minutes)) + delim + str(seconds)
 
# function to convert degrees into dec deg:mm:ss string
def _deg_to_dec_dms(angle_deg, delim = ':', round_seconds = False):
    if angle_deg >= 0:
        deg = angle_deg - (angle_deg % 1)
        minutes = (angle_deg % 1) * 60
    else:
        deg = angle_deg + (-angle_deg % 1)
        minutes = (-angle_deg % 1) * 60
 
    seconds = (minutes % 1) * 60
     
    # decide to round the seconds (for display) or leave unrounded (internal)
    if round_seconds: return str(int(deg)) + delim + str(int(minutes)) + delim + str(round(seconds))
    else: return str(int(deg)) + delim + str(int(minutes)) + delim + str(seconds)
 
# function to get image urls from ESO database given ra dec in degrees
def _ra_dec_img_url(ra_deg, dec_deg, imageonly = True):
    # convert to hms/dms
    ra = _deg_to_ra_hms(ra_deg, delim = '%3A')
    dec = _deg_to_dec_dms(dec_deg, delim = '%3A')
     
    # find first url and read off html
    url_full = 'https://archive.eso.org/dss/dss/image?'
    url_full += 'ra=' + ra
    url_full += '&dec=' + dec
    url_full += '&equinox=J2000&name=&x=5&y=5&Sky-Survey=DSS1&mime-type=image%2Fgif&statsmode=WEBFORM'
     
    # search for true image url
    if imageonly:
        url_full_HTML = str( urllib.request.urlopen(url_full).read() )
        
        start = url_full_HTML.find('IMG SRC=/dss/dss') + 8
        end = url_full_HTML.find('.gif ALT') + 4
        url_pic = 'https://archive.eso.org'
        for i in range(start, end):
            url_pic += url_full_HTML[i]
         
        return url_pic
     
    # want full url, not just the image (which expires anyway)
    else: return url_full
    
# function to write the html page
def htmlwrite(exoplanet_array, tele_array, filename, print_properties, path):
    # output to html file
    outfile = open(filename + '.html','w')
    outfile.write("<html>\n<head>\n")
     
    # table style
    style = "<style>"
    style += """
    table {
        border-collapse: collapse; /* collapsed borders */
        width: 100%; /* full width */
        }\n"""
    style += "table, th, td { border: 1px solid black; } /* black borders */\n" 
    style += """
    th {
        background-color: #4CAF50; /* green background for header */
        color: white; /* white header text */
        padding: 12px; /* header padding */
        font-size: 18px; /* font size */
        }\n"""
    style += "tr:nth-child(even) { background-color: #f2f2f2; }\n /* every other row shaded grey */"
    style += ".fixed_header thead th { position: sticky; top: 0; }\n /* frozen header */"
     
    # search button styling by https://www.w3schools.com/howto/howto_js_filter_table.asp
    style += """
    #userInput {
      background-repeat: no-repeat; /* Do not repeat the icon image */
      width: 100%; /* Full-width */
      font-size: 16px; /* Increase font-size */
      padding: 12px 20px 12px 40px; /* Add some padding */
      border: 1px solid #ddd; /* Add a grey border */
      margin-bottom: 12px; /* Add some space below the input */
      }\n"""
    style += "</style>\n"
    outfile.write(style)
     
    outfile.write("</head>\n<body>")
    outfile.write("""<input type="text" id="userInput" onkeyup="searchFunction()" placeholder="Search for names..." title="Type in a name">""")
     
    # create table
    table = """<table class="fixed_header" id="planetTable">\n"""
     
    # table header
    table += "<thead>\n<tr>\n"
    table += "<th>image</th>\n"
    for key in print_properties:
        table += "<th>" + key.replace('_', ' ') + "</th>\n"
    table += "</tr>\n</thead>\n"
     
    # table body
    table += "<tbody>\n"
    for planet in exoplanet_array:
        table += "<tr>\n"
        
        # find image by url if necessary but store the image for indefinite use
        foundImage = False
        url_full = _ra_dec_img_url(planet.ra, planet.dec, imageonly = False)
         
        # search for image file in directory
        for imagefilename in os.listdir(path):
            if imagefilename == planet.name.replace(' ', '_').replace('/', '_') + '.gif':
                # image found in file, put in table
                foundImage = True
                table += """<td><a href="%s"><img src="%s%s" alt = "%s" style="border:0"></a></td>\n""" %(url_full, path, imagefilename, planet.name)
                 
        # if no image found, download a new one
        if not foundImage:
            # download new image
            url_pic = _ra_dec_img_url(planet.ra, planet.dec)
            imagefilename = path + planet.name.replace(' ', '_').replace('/', '_') + '.gif'
            urllib.request.urlretrieve(url_pic, imagefilename)
             
            # file downloaded, put in table
            table += """<td><a href="%s"><img src="%s%s" alt = "%s" style="border:0"></a></td>\n""" %(url_full, path, imagefilename, planet.name)
        
        # print all other properties
        for key in print_properties:
            # if printing the name, make it a hyperlink
            if key == 'name':
                link = "http://www.exoplanet.eu/catalog/%s/" %planet.name.replace(' ', '_')
                display = planet.name
                table += "<td><a href=\"%s\">%s" %(link, display) + "</td>\n"
             
            # for ra and dec display in hms/dms
            elif key == 'ra':
                table += "<td>" + _deg_to_ra_hms(planet.ra, round_seconds = True) + "</td>\n"
            elif key == 'dec':
                table += "<td>" + _deg_to_dec_dms(planet.dec, round_seconds = True) + "</td>\n"
             
            # otherwise display as normal
            else:
                table += "<td>" + str(getattr(planet, key)) + "</td>\n"
         
        table += "</tr>\n"
     
    table += "</tbody>\n</table>"
     
    outfile.write(table)
     
    # search function from https://www.w3schools.com/howto/howto_js_filter_table.asp
    outfile.write("""
    <script>
    function searchFunction() {
      // Declare variables 
      var input, filter, table, tr, td, i, txtValue;
      input = document.getElementById("userInput");
      filter = input.value.toUpperCase();
      table = document.getElementById("planetTable");
      tr = table.getElementsByTagName("tr");
     
      // Loop through all table rows, and hide those who don't match the search query
      for (i = 0; i < tr.length; i++) {
        td = tr[i].getElementsByTagName("td")[1]; // element 1 for name column
        if (td) {
          txtValue = td.textContent || td.innerText;
          if (txtValue.toUpperCase().indexOf(filter) > -1) {
            tr[i].style.display = "";
          } else {
            tr[i].style.display = "none";
          }
        } 
      }
    }
    </script>\n""")
     
    outfile.write("</body>\n</html>")
     
    outfile.close()