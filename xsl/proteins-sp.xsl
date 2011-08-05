<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>


    <xsl:template match="proteins">
        <html>
<head>
            <title>Signal peptides search page</title>
    <script type="text/javascript" src="acids.js"/>
    <style>
        span.match,td.break{font-weight:bold;};
        {font-weight:bold;};
    </style>

</head>
            <body>
                <script>
                    function getMass(m, mod) {
                        if (mod == "ammonia") {
                            m += ammonia;
                        }
                        if (mod == "aw") {
                            m += ammonia + water;
                        }
                        if (mod == "ww") {
                            m += water + water;
                        }
                        if (mod == "water") {
                            m += water;
                        }
                        if (mod == "plus") {
                            m += 1;
                        }
                        if (mod == "minus") {
                            m -= 1;
                        }
                        return m;
                    }


                    function acidMass(a) {
                        for (var i = 0; acids.length > i; i++) {
                                var acid = acids[i];
                                if (acid.letter == a) {
                                    return acid.monoMass;
                                }
                        }
                    }

                    var min = -999999999;

                    var modes = [null, "plus", "minus",  "ammonia", "water",  "aw", "ww"]

                    var proton = 1.007276;
                    var isotope = 1.00235;
                    var water = 18.010565;
                    var ammonia = 17.0265;
                    var oxygen = 15.9949;
                    var CO = oxygen + 12;
                    var H = 1.007825;
                    var N = 14.0030740048;
                    var NH = N + H;

                    function check(p, total,  t) {
                        var r = total - p - water;
                        var m = [
                                p,
                                r,
                                p-1, p+1, p + CO,
                                r-1, r+1, r - ammonia, r - water, r + water
                            ];
                            //24,0,1,13,1,1,4,6,12,8,3,3
                        for(var i = 0; m.length > i; i++) {
                            if (0.1 > Math.abs(m[i] - t)) {
                                stat[i]++;
                                return m[i];
                            }
                        }
                        return min - 1;
                    }

                    var stat = [];


                    function renderMatch(protein, spectrum, shift, mass) {
                        stat = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
                        var text = "&lt;font color='grey'>";
                        var targetMass = -shift ;
                        var grey = true;
                        for (var cur = 0; protein.length > cur; cur++) {
                            var add = "";
                            if (grey &amp;&amp; cur%80 == 0 ) {
                                if (cur > 0) {
                                    add += "&lt;br/>";
                                }
                            }
                            for (var i = 0;  spectrum.length > i; i++) {
                                var v = spectrum[i];
                                var vopt = check(v, mass, targetMass);
                                if (vopt > min) {
                                    if (grey)  {
                                        grey = false;
                                        add += "&lt;/font>&lt;br/>&lt;font color='red'>" + vopt  + "&lt;/font>&lt;br/>";
                                    }
                                    add += "|";
                                }
                            }
                            text += add + protein.charAt(cur);
                            targetMass += acidMass(protein.charAt(cur));
                        }
                        //text = stat + text
                        return text;
                    }

                </script>
                <xsl:apply-templates select="protein">
                    <xsl:sort select="protein-id" data-type="number"/>
                </xsl:apply-templates>

            </body>
        </html>
    </xsl:template>

    <xsl:template match="protein">
        <h4>Protein #<xsl:value-of select="protein-id"/> - <xsl:value-of select="name"/></h4>

        <script>
            var protein<xsl:value-of select="protein-id"/> = '<xsl:apply-templates select="acids"/>';
        </script>

        <xsl:apply-templates select="matches/match[first-residue > 0]/spectrum">
            <xsl:sort select="../n-mass" data-type="number"/>
        </xsl:apply-templates>

    </xsl:template>


    <xsl:template match="spectrum">
        <div>Match with spectrum #<a href="http://msalign.ucsd.edu/set8-7/html/prsms/prsm{spectrum-id}.html"><xsl:value-of select="spectrum-id"/></a>
            with score  <xsl:value-of select="../best-score"/></div>
        <div id="spectrum{spectrum-id}">
        </div>
        <!--<div id="spectrum{spectrum-id}"></div>-->
        <script>
            var spectrum<xsl:value-of select="spectrum-id"/> =  [ 0
        <xsl:apply-templates select="peaks/peak"/>
            ];

            var s =
            "&lt;font color='grey'>" + protein<xsl:value-of select="../../../protein-id"/>.substr(0, <xsl:value-of select="../first-residue"/>) +
            "&lt;/font>&lt;br/>&lt;font color='red'>" + <xsl:value-of select="../n-mass"/>  + "&lt;/font>&lt;br/>"+
            protein<xsl:value-of select="../../../protein-id"/>.substr(<xsl:value-of select="../first-residue"/>);
            //document.getElementById("spectrumj<xsl:value-of select="spectrum-id"/>").innerHTML = s;
            document.getElementById("spectrum<xsl:value-of select="spectrum-id"/>").innerHTML = renderMatch(
                protein<xsl:value-of select="../../../protein-id"/>,
                spectrum<xsl:value-of select="spectrum-id"/>,
                <xsl:value-of select="../best-shift"/>,
                <xsl:value-of select="precursor-mass"/>
            );
        </script>
    </xsl:template>

    <xsl:template match="peak">
        , <xsl:value-of select="monoisotopic-mass"/>
    </xsl:template>

    <xsl:template match="residue">
        <xsl:param name="break"/>
        <xsl:if test="residue-id = $break or residue-id = $break + 1">
            <xsl:value-of select="acid"/>
        </xsl:if>
    </xsl:template>

    <xsl:template match="residue" mode="full">
        <xsl:param name="break"/>
        <xsl:if test="(residue-id > ($break -3)) and  (($break + 4)>residue-id)">
            <xsl:value-of select="acid"/>
        </xsl:if>
    </xsl:template>

    <xsl:template match="acid"><xsl:value-of select="."/></xsl:template>
</xsl:stylesheet>
