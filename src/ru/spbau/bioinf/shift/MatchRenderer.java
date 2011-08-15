package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import org.jdom.Document;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.util.List;
import java.util.Map;

public class MatchRenderer {

    private static final Logger log = Logger.getLogger(MatchRenderer.class);

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);
        MatchRenderer renderer = new MatchRenderer();
        renderer.render(config);
    }

    public void render(Configuration config) throws Exception {
        log.debug("star result processing");
        Map<Integer, Spectrum> spectrums = config.getSpectrums();

        log.debug("spectrums data loaded");
        List<Protein> proteins = config.getProteins();

        BufferedReader matchReader = ReaderUtil.getBufferedReader(config.getMatchFile());
        String s;
        ScoringFunction scoringFunction = new ExtendedSharedPeaksScoringFunction();
        while ((s = matchReader.readLine()) != null) {
            String[] data = s.split(" ");
            int spectrumId = Integer.parseInt(data[0]);
            int proteinId = Integer.parseInt(data[1]);
            Spectrum spectrum = spectrums.get(spectrumId);
            Protein protein = proteins.get(proteinId);
            SpectrumProteinMatch spm = new SpectrumProteinMatch(spectrum, protein, scoringFunction);
            Document doc = new Document();
            doc.setRootElement(spm.toXml());
            XmlUtil.saveXml(doc, config.getSpectrumXmlFile(spectrum));
        }
    }

}
