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

                    var proton = 1.007276;
                    var isotope = 1.00235;
                    var water = 18.010565;
                    var ammonia = 17.0265;
                    var oxygen = 15.9949;


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

                    function check(v, total,  t) {
                        for(var i = 0; modes.length > i; i++) {
                            var m = getMass(v, modes[i]);
                            if (0.1 > Math.abs(m - t)) {
                                return m;
                            }
                            if (0.1 > Math.abs(total - m - t)) {
                                return total - m;
                            }
                        }
                        return min - 1;
                    }


                    function renderMatch(protein, spectrum, shift, mass) {
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
        <h2>Protein #<xsl:value-of select="protein-id"/> match</h2>

        <script>
            var protein<xsl:value-of select="protein-id"/> = '<xsl:apply-templates select="acids"/>';
        </script>

        <xsl:apply-templates select="matches/match/spectrum">
            <xsl:sort select="../first-residue" data-type="number"/>
        </xsl:apply-templates>

    </xsl:template>


    <xsl:template match="spectrum">
        <div>Match with spectrum #<a href="spectrums/spectrum{spectrum-id}.html"><xsl:value-of select="spectrum-id"/></a></div>
        <div id="spectrumj{spectrum-id}">
        </div>
        <div id="spectrum{spectrum-id}"></div>
        <script>
            var spectrum<xsl:value-of select="spectrum-id"/> =  [ 0
        <xsl:apply-templates select="peaks/peak"/>
            ];
            var s =
            "&lt;font color='grey'>" + protein<xsl:value-of select="../../../protein-id"/>.substr(0, <xsl:value-of select="../first-residue"/>) +
            "&lt;/font>&lt;br/>&lt;font color='red'>" + <xsl:value-of select="../n-mass"/>  + "&lt;/font>&lt;br/>"+
            protein<xsl:value-of select="../../../protein-id"/>.substr(<xsl:value-of select="../first-residue"/>);
            document.getElementById("spectrumj<xsl:value-of select="spectrum-id"/>").innerHTML = s;

        </script>
        <a href="#" onclick="document.getElementById('spectrum{spectrum-id}').innerHTML=renderMatch(protein{../../../protein-id}, spectrum{spectrum-id},  {../best-shift}, {precursor-mass}); return false;">show breaks</a>
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
