package ru.spbau.bioinf.shift;

public interface ScoringFunction {
    double getScore(Spectrum spectrum, Protein protein, double shift);
}
