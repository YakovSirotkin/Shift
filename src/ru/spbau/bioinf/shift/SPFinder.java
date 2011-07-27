package ru.spbau.bioinf.shift;

import org.jdom.Document;
import org.jdom.Element;
import ru.spbau.bioinf.shift.util.ReaderUtil;
import ru.spbau.bioinf.shift.util.XmlUtil;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SPFinder {

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);

        ProteinFinder finder = new ProteinFinder(config);
        List<Spectrum> spectrums = finder.getSpectrums();
        List<Protein> proteins = finder.getProteins();

        BufferedReader matchReader = ReaderUtil.getBufferedReader(config.getMatchFile());

        HashMap<Integer, List<Integer>> hits = new HashMap<Integer, List<Integer>>();
        String s;
        while ((s = matchReader.readLine()) != null) {
            String[] data = s.split(" ");
            int spectrumId = Integer.parseInt(data[0]);
            int proteinId = Integer.parseInt(data[1]);
            if (!hits.containsKey(proteinId)) {
                hits.put(proteinId, new ArrayList<Integer>());
            }
            hits.get(proteinId).add(spectrumId);
        }

        Document doc = new Document();
        Element root = new Element("proteins");
        doc.setRootElement(root);

        for (Map.Entry<Integer, List<Integer>> entry : hits.entrySet()) {
            Protein protein = proteins.get(entry.getKey());
            List<Integer> matches = entry.getValue();
            for (int spectrumId : matches) {
                protein.addSpectrumMatch(spectrums.get(spectrumId));
            }
            root.addContent(protein.toXml());
        }

        XmlUtil.saveXml(doc, config.getProteinXmlFile());
    }
}
