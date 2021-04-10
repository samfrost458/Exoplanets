function deg_to_ra_hms(angle_deg, delim = ':', round_seconds = false) {
    var hours, minutes, seconds;
    hours = (angle_deg - (angle_deg % 15)) / 15
    minutes = (angle_deg / 15 - hours) * 60
    seconds = (minutes % 1) * 60
     
    // decide to round the seconds (for display) or leave unrounded (internal)
    if (round_seconds) {
        return String(Math.floor(hours)) + delim + String(Math.floor(minutes)) + delim + String(seconds.toFixed());
    }
    else {
        return String(Math.floor(hours)) + delim + String(Math.floor(minutes)) + delim + String(seconds);
    }
}

function deg_to_dec_dms(angle_deg, delim = ':', round_seconds = false) { 
    var deg, minutes, seconds;   
    if (angle_deg >= 0) {
        deg = angle_deg - (angle_deg % 1);
        minutes = (angle_deg % 1) * 60;
    }
    else {
        deg = angle_deg + (-angle_deg % 1);
        minutes = (-angle_deg % 1) * 60;
    }
 
    seconds = (minutes % 1) * 60;
     
    // decide to round the seconds (for display) or leave unrounded (internal)
    if (round_seconds) {
        return String(Math.floor(deg)) + delim + String(Math.floor(minutes)) + delim + String(seconds.toFixed());
    }
    else {
        return String(Math.floor(deg)) + delim + String(Math.floor(minutes)) + delim + String(seconds);
    }
}

function ra_dec_img_url(ra_deg, dec_deg) {
    var ra, dec, url_full;
    // convert to hms/dms
    ra = deg_to_ra_hms(ra_deg, delim = '%3A');
    dec = deg_to_dec_dms(dec_deg, delim = '%3A');
     
    // find first url and read off html
    url_full  = 'https://archive.eso.org/dss/dss/image?';
    url_full += 'ra=' + ra;
    url_full += '&dec=' + dec;
    url_full += '&equinox=J2000&name=&x=5&y=5&Sky-Survey=DSS1&mime-type=image%2Fgif&statsmode=WEBFORM';
     
    return url_full;
}

function round_sig_figs(number, significant) {
    // try for float
    try { return number.toFixed(significant); }
    // if it doesn't work number is probably a string
    catch(err) { return number; }
}

function capitalize(str) {
  return str.charAt(0).toUpperCase() + str.slice(1)
}

