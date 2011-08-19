<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>


    <xsl:template match="match">
        <html>
            <xsl:apply-templates select="spectrum"/>
        </html>
    </xsl:template>

    <xsl:template match="spectrum">
            <head>
                <title>Spectrum #<xsl:value-of select="spectrum-id"/> match</title>
                <script type="text/javascript" src="../acids.js"/>
                <style>
                    span.match,td.break{font-weight:bold;}
                    a.peak{ text-decoration: none;}
                </style>
            </head>
            <body>

                <xsl:call-template name="navigation"/>

                <h2>Spectrum #<xsl:value-of select="spectrum-id"/> match</h2>

                <table cellpadding="3">
                    <tr>
                        <td>
                             Precursor m/z:
                        </td>
                        <td>
                            <xsl:value-of select="precursor-mz"/>
                        </td>
                        <td>Precursor charge:</td>
                        <td>
                            <xsl:value-of select="precursor-charge"/>
                        </td>
                        <td>Precursor mass:</td>
                        <td>
                            <xsl:value-of select="precursor-mass"/>
                        </td>
                    </tr>
                    <tr>
                        <td colspan="8">
                            Scan(s): <xsl:value-of select="scans"/>
                        </td>
                    </tr>
                </table>
                <br/>

                <script>
                    var proton = 1.007276;
                    var isotope = 1.00235;
                    var water = 18.010565;
                    var ammonia = 17.0265;
                    var oxygen = 15.9949;
                    var CO = oxygen + 12;
                    var H = 1.007825;
                    var N = 14.0030740048;
                    var NH = N + H;

                    var proteins = [];

                    var total = <xsl:value-of select="precursor-mass"/> - 0.035;
                    var peaks = [
                                {
                                    id: -1,
                                    mass: 0,
                                    intensity: 0,
                                    attach:new Array(),
                                    selected: [],
                                    reverse: 0
                                },
                        <xsl:apply-templates select="peaks/peak"/>
                                {
                                    id: -2,
                                    mass: total,
                                    intensity: 0,
                                    attach:new Array(),
                                    selected: [],
                                    reverse: 0
                                },
                    ];

                    var all = [];

                    for (var i = 0; peaks.length > i; i++) {
                        all[peaks[i].id] = peaks[i];
                    }

                </script>
                <script>
                    function acidMass(a) {
                        for (var i = 0; acids.length > i; i++) {
                                var acid = acids[i];
                                if (acid.letter == a) {
                                    return acid.monoMass;
                                }
                        }
                    }

                     var br = "&lt;br&gt;\r";

                    function compareScores(s1, s2) {
                        return s2.score - s1.score;
                    }

                    function unclear(p, level) {
                        for (var i = 0; level >= i; i++) {
                            if (p.selected[i]) {
                                return false;
                            }
                        }
                        return true;
                    }

                    var protein;

                    var sharedPeaks;
                    var describedPeaks;

                    function renderMatch(proteinId, shift) {
                        protein = proteins[proteinId];
                        sharedPeaks = 0;
                        describedPeaks = 0;
                        var text = "";
                        var targetMass = -shift ;
                        for (var i = 0; peaks.length > i; i++) {
                            var p = peaks[i];
                            p.dist = 100000000;
                            p.place = -1;
                            var targetMass = -shift;
                            for (var cur = 0; protein.length > cur; cur++) {
                                if (Math.abs(p.dist) > Math.abs(p.mass - targetMass)) {
                                    p.dist =  p.mass - targetMass;
                                    p.place = cur;
                                    p.reverse = 0;
                                }
                                if (p.id > 0) {
                                    var r =  total - p.mass;
                                    if (Math.abs(p.dist) > Math.abs(r - targetMass)) {
                                        p.dist =  r - targetMass;
                                        p.place = cur;
                                        p.reverse = 1;
                                    }
                                }
                                targetMass += acidMass(protein.charAt(cur));
                            }
                        }
                        var targetMass = -shift ;
                        unmatchedPeaks = [];
                        for (var cur = 0; protein.length > cur; cur++) {
                            var b = [];
                            for (var i = 0; peaks.length > i; i++) {
                                var p = peaks[i];
                                if (p.place == cur) {
                                    b[b.length] = p;
                                }
                            }
                            if (b.length > 0) {
                                text += br + "&lt;span class='peak'>";
                                b.sort(compareBreaks);
                                for (var i = 0; b.length > i; i++) {
                                    text += "&amp;#160;";
                                    var title = b[i].mass;
                                    if (b[i].reverse==1) {
                                        title = "(" + (total - b[i].mass)  + ") " + title;
                                    }
                                    title+=" " +  b[i].dist;
                                    var link = "&lt;a href='#' class='peak' title='" + title + "'>" + getText(b[i]) + "&lt;/a>";
                                    text += link;
                                }
                                text += "&lt;/span>";
                                text += br;
                            }
                            text += protein.charAt(cur);
                            targetMass += acidMass(protein.charAt(cur));
                        }
                        text +=  br + br;
                        document.getElementById("p" + proteinId).innerHTML = text;
                        document.getElementById("message").innerHTML =
                            sharedPeaks + " matched peaks,  " + describedPeaks + " described peaks, " + unmatchedPeaks.length+ " unmatched peaks.";

                        var table = document.getElementById("unmatched");
                        while (table.rows.length > 1) {
                            table.deleteRow(1);
                        }

                        unmatchedPeaks.sort(compareMass);
                        for (var i = 0; unmatchedPeaks.length > i; i++) {
                            var p = unmatchedPeaks[i];
                            var row = table.insertRow(i+1);
                            var c = [];
                            for (var j = 0; 13 > j; j++) {
                                c[j] = row.insertCell(j);
                            }

                            c[0].innerHTML = p.id;
                            renderMass(p.mass, c, 1, shift);
                            if (p.id > 0) {
                                renderMass(total - p.mass, c, 7, shift);
                            }
                        }

                        var g = [];
                    for (var i = 0; unmatchedPeaks.length > i; i++){
                        var p = unmatchedPeaks[i];
                            g[g.length] = {mass: p.mass, reverse: 0, id: p.id};
                            if (p.id > 0) {
                                g[g.length] = {mass: total - p.mass, reverse: 1, id: p.id};
                            }
                        }
                        g.sort(compareMass);
                        table = document.getElementById("graph");
                        while (table.rows.length > 1) {
                            table.deleteRow(1);
                        }
                        var prev = -100000;
                        for (var i = 0; g.length > i; i++) {
                            var v = g[i];
                            var row = table.insertRow(i+1);
                            var c = [];
                            for (var j = 0; 2 > j; j++) {
                                c[j] = row.insertCell(j);
                            }
                            c[0].innerHTML = v.id;
                            c[1].innerHTML = v.reverse == 0 ? v.mass : "(" + v.mass + ")";
                            var pos = 2;
                            for (var j = i; g.length > j; j++) {
                                var diff = getAcid(g[j].mass - prev);
                                if (diff. length == 1) {
                                    c[pos] = row.insertCell(pos);
                                    c[pos].rowSpan = j + 1 - i;
                                    c[pos].innerHTML = diff;
                                    pos++;
                                }
                            }

                            prev = v.mass;
                        }
                    }

                    function renderMass(m, c, start, shift) {
                        c[start].innerHTML = m.toFixed(3);
                        var targetMass = -shift ;
                        for (var cur = 0; protein.length > cur; cur++) {
                            var d = acidMass(protein.charAt(cur));
                            if (targetMass + d > m) {
                                c[start + 2].innerHTML = getAcid(m - targetMass);
                                var s = "";
                                s += cur > 0 ? protein.charAt(cur - 1) : "_";
                                if (cur > 0) {
                                    c[start + 1].innerHTML = getAcid(m - targetMass + acidMass(protein.charAt(cur - 1)));
                                }
                                s += protein.charAt(cur);
                                s += protein.length - 1 > cur ? protein.charAt(cur + 1) : "_";
                                if (protein.length - 1 > cur) {
                                    c[start +5].innerHTML = getAcid(- m + targetMass + d + acidMass(protein.charAt(cur + 1)));
                                }
                                c[start + 3].innerHTML = s;
                                c[start + 4].innerHTML = getAcid(- m + d + targetMass) ;
                                break;
                            }
                            targetMass += d;
                        }
                    }

                    function getAcid(v) {
                        for (var i = 0; acids.length > i; i++) {
                            var acid = acids[i];
                            if (0.02 > Math.abs(acid.monoMass -v)) {
                                return acid.letter;
                            }
                        }
                        return v.toFixed(3);
                    }

                    var unmatchedPeaks = [];

                    function compareBreaks(p1, p2) {
                        return p1.dist - p2.dist;
                    }

                    function compareMass(p1, p2) {
                        return p1.mass - p2.mass;
                    }

                    function processShifts(proteinId) {
                        for (var i =  0; shifts.length > i; i++) {
                            var shift = shifts[i];
                        }
                        shifts.sort(compareScores);
                        var d = document.getElementById("shifts" + proteinId);
                        text = "";
                        for (var i =  0; shifts.length > i; i++) {
                            var shift = shifts[i];
                            if (shift.score > 0) {
                                text +='<a href="#" onclick="renderMatch(' + shift.proteinId + ', '  +  shift.value + ');">' + shift.value.toFixed(3) + '(' + shift.score+')</a> ';
                            }
                        }
                        d.innerHTML = text;
                        renderMatch(shifts[0].proteinId, shifts[0].value)
                    }

                    var bmod = [
                    {value: -1, text: "-1"},
                    {value: 1, text: "+1"},
                    {value: -CO, text: "-CO"},
                    {value: -water, text: "-water"},
                    {value: water, text: "+water"},
                    ];
                    var ymod = [
                        {value: -1, text: "(+1)"},
                        {value: 1, text: "(+1)"},
                        {value: ammonia, text: "(-ammonia)"},
                        {value: -water, text: "(+water)"},
                        {value: water, text: "(-water)"},
                    ];

                    var epsilon = 0.02;

                    function getText(p) {
                        var color="red";
                        var text= p.dist.toFixed(2);
                        if (p.reverse==1) {
                            text = "(" + text + ")";
                        }
                        if (epsilon > Math.abs(p.dist)) {
                            color="green";
                            text = p.reverse == 1 ? "y" : "b";
                            if (p.id == -1) {
                                text = "begin";
                            }
                            if (p.id == -2) {
                                text = "end";
                            }
                        }
                        if (p.id > 0) {
                            if (p.reverse==0) {
                                for (var i =0; bmod.length > i; i++) {
                                    if (epsilon > Math.abs(p.dist - bmod[i].value)) {
                                        color = "yellowgreen";
                                        text = bmod[i].text;
                                    }
                                }
                            }
                            if (p.reverse==1) {
                                for (var i =0; ymod.length > i; i++) {
                                    if (epsilon > Math.abs(p.dist - ymod[i].value)) {
                                        color = "yellowgreen";
                                        text = ymod[i].text;
                                    }
                                }
                            }
                        }

                        if (color == 'red') {
                            unmatchedPeaks[unmatchedPeaks.length] = p;
                        } else if (color == 'green') {
                            sharedPeaks ++;
                        } else {
                            describedPeaks ++;
                        }

                        return "&lt;font color='" + color + "'>" + text + "&lt;/font>";
                    }

                </script>


                <h3>Match with protein #<xsl:value-of select="../protein/protein-id"/> &#160;  <xsl:value-of select="../protein/name"/></h3>

                <h3>Shifts</h3>
                <div id="shifts{../protein/protein-id}">
                </div>

                <h3>Protein with peaks associated with the nearest cleavages</h3>
                <p id="message"></p>
                <p id="p{../protein/protein-id}" style="font-family: 'FreeMonoMedium', FreeMonoMedium, Miltonian, monospace;"></p>

                <h3>Distances between unmatched  peaks and the nearest cleavages in protein</h3>
                <table id="unmatched" border="1">
                    <tr><th>Peak Id</th><th>b-ion</th><th colspan="5">b position</th><th>y-ion</th><th colspan="5">y-position</th></tr>
                </table>
                <h3>Amino acids between unmatched peaks</h3>
                <table id="graph" border="1">
                    <tr><th>Peak Id</th><th>Mass</th></tr>
                </table>
                <script>
                    proteins[<xsl:value-of select="../protein/protein-id"/>] = '<xsl:value-of select="../protein/acids"/>';
                    var shifts =  [
                <xsl:apply-templates select="../shifts/shift"/>
                    ];
                    processShifts(<xsl:value-of select="../protein/protein-id"/>,0);
                </script>

                <xsl:call-template name="navigation"/>
            </body>
    </xsl:template>

    <xsl:template match="shift">
        {proteinId:<xsl:value-of select="../../protein/protein-id"/>, value:<xsl:value-of select="value"/>, score:<xsl:value-of select="score"/>}
        <xsl:if test="last() > position()">,</xsl:if>
    </xsl:template>

    <xsl:template match="peak">
        {
            id: <xsl:value-of select="position()"/>,
            mass: <xsl:value-of select="monoisotopic-mass"/>,
            intensity: <xsl:value-of select="intensity"/>,
            attach:new Array(),
            selected: [],
            reverse: 0
        },
    </xsl:template>

    <xsl:template name="navigation">
        <p>
             <a href="../proteins-sp.html">Signal peptides</a>
        </p>
    </xsl:template>
</xsl:stylesheet>
