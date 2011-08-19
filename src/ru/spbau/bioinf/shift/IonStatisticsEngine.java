package ru.spbau.bioinf.shift;

import org.apache.log4j.Logger;
import ru.spbau.bioinf.shift.util.ReaderUtil;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class IonStatisticsEngine {

    private static final Logger log = Logger.getLogger(IonStatisticsEngine.class);

    public static void main(String[] args) throws Exception {
        Configuration config = new Configuration(args);
        IonStatisticsEngine stat = new IonStatisticsEngine();
        stat.render(config);
    }

    int limit = 100 * 100;
    int[] statY = new int[limit];

    public static double getBestDiff(List<Double> sharedPeaks, double m) {
        double bestDiff = 10000;
        for (double sharedPeak : sharedPeaks) {
            if (Math.abs(bestDiff) > Math.abs(m - sharedPeak)) {
                bestDiff = m - sharedPeak;
            }
        }
        return bestDiff;
    }

    public void render(Configuration config) throws Exception {
        log.debug("star result processing");
        Map<Integer, Spectrum> spectrums = config.getSpectrums();

        log.debug("spectrums data loaded");
        List<Protein> proteins = config.getProteins();

        BufferedReader matchReader = ReaderUtil.getBufferedReader(new File(config.getMatchFile().getParent(), "match_bonly.txt"));
        String s;
        List<double[]> points = new ArrayList<double[]>();
        while ((s = matchReader.readLine()) != null) {
            String[] data = s.split(" ");
            int spectrumId = Integer.parseInt(data[0]);
            int proteinId = Integer.parseInt(data[1]);
            double shift = Double.parseDouble(data[3]);
            Spectrum spectrum = spectrums.get(spectrumId);
            Protein protein = proteins.get(proteinId);
            double[] pd = protein.getSpectrum();
            List<Double> sharedPeaks = new ArrayList<Double>();
            for (double v : pd) {
                sharedPeaks.add(v - shift);
            }
            spectrum.clearData();
            List<Peak> peaks = spectrum.getPeaks();
            for (Peak peak : peaks) {
                double m = peak.getMonoisotopicMass();
                    double bestDiff = getBestDiff(sharedPeaks, spectrum.getPrecursorMass() - m);
                    int val = (int) Math.round(bestDiff * 1000);
                    int pos = val;// + 100;
                    if (pos >= 0 && pos < limit) {
                        statY[pos]++;
                    }
                    if (bestDiff > -0.02 && bestDiff <= 0.05) {
                        points.add(new double[]{m, bestDiff});
                    }
            }
        }
        double sum = 0;
        for (double[] point : points) {
            sum+= point[1];
        }

        System.out.println("Stat for Y-ions");
        for (int i = 0; i < limit; i++) {
            System.out.println(i + " " + statY[i]);
        }

        findApproximation(points, sum);
    }

    private void findApproximation(List<double[]> points, double sum) {
        double average = //0.03655;
        sum / points.size();
        System.out.println("average = " + average);
        System.out.println("points.size() = " + points.size());

        double minA = 0.02;
        double maxA = 0.05;
        double maxK = 0.000001;
        double bestDiff = getDiff(points, average, 0);
        System.out.println("bestDiff = " + bestDiff);
        Random r = new Random();
        System.out.println("maxK = " + maxK);
        while(2>maxK) {
            double da = r.nextDouble();
            double a = minA + da * (maxA - minA);
            double dk = r.nextDouble();
            double k = -maxK + dk * 2 * maxK;
            double nextDiff = getDiff(points, a, 0);
            if (nextDiff < bestDiff) {
                bestDiff = nextDiff;
                System.out.println(a + " " + " " + k + " " + nextDiff);
            }
        }
    }

    private double getDiff(List<double[]> points, double average, double k) {
        double diff = 0;
        for (double[] point : points) {
            diff += Math.abs(point[1] - average - k * point[0]);
        }
        return diff;
    }
}
