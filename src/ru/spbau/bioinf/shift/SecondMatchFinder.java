package ru.spbau.bioinf.shift;

import ru.spbau.bioinf.shift.util.ReaderUtil;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class SecondMatchFinder {

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);
        ScoringFunction scoringFunction = config.getScoringFunction();
        ProteinFinder finder = new ProteinFinder(config);
        List<Spectrum> spectrums = finder.getSpectrums();
        List<Protein> proteins = finder.getProteins();

        BufferedReader matchReader = ReaderUtil.getBufferedReader(config.getMatchFile());

        PrintWriter matchFile = ReaderUtil.createOutputFile(config.getSecondMatchFile());
        PrintWriter matchPairsFile = ReaderUtil.createOutputFile(config.getPairsMatchFile());

        String s;
        while ((s = matchReader.readLine()) != null) {
            String[] data = s.split(" ");
            int spectrumId = Integer.parseInt(data[0]);
            int proteinId = Integer.parseInt(data[1]);
            Spectrum sp = spectrums.get(spectrumId);
            Protein p = proteins.get(proteinId);
            double[] pd = p.getSpectrum();
            Map<String,List<Double>> positions = ProteinFinder.getPositions(sp);
            List<Double> shifts = ProteinFinder.getShifts(p, positions);
            double maxScore = 0;
            double bestShift = 0;
            for (double shift : shifts) {
                double score = scoringFunction.getScore(sp, p, shift);
                if (score > maxScore) {
                    maxScore = score;
                    bestShift = shift;
                }
            }
            Spectrum left = sp.getLeft(pd, bestShift);
            int totalPeaksCount = sp.getPeaks().size();
            maxScore = totalPeaksCount - left.getPeaks().size();
            positions = ProteinFinder.getPositions(left);
            shifts = ProteinFinder.getShifts(p, positions);
            double maxScoreLeft = 0;
            double bestShiftLeft = 0;
            for (double shift : shifts) {
                double score = scoringFunction.getScore(left, p, shift);
                if (score > maxScoreLeft) {
                    maxScoreLeft = score;
                    bestShiftLeft = shift;
                }
            }
            Spectrum leftLeft = left.getLeft(pd, bestShiftLeft);
            maxScoreLeft = left.getPeaks().size() - leftLeft.getPeaks().size();
            matchFile.println(spectrumId + " " + proteinId + " " + maxScore + " " + bestShift + " " + maxScoreLeft + " " +  bestShiftLeft + " " + totalPeaksCount);
            matchFile.flush();
            ArrayList<Match> matches = finder.getSpectrumMatches(left);
            if (matches.size() == 1) {
                Match match = matches.get(0);
                matchPairsFile.println(spectrumId + " " +  proteinId + " " + maxScore + " " + match.getProteinId() + " " + match.getScore() + " " + totalPeaksCount);
                matchPairsFile.flush();
            }


        }
        matchFile.close();
        matchPairsFile.close();
    }
}
