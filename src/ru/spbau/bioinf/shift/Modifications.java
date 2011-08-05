package ru.spbau.bioinf.shift;

import java.util.Arrays;
import java.util.List;

public enum Modifications {

    CORE(){
        @Override
        List<Double> getModifications(double p, double precursorMass) {
            double r = precursorMass - p - Consts.WATER;
            return Arrays.asList(p, r);
        }
    },
    ADDITIONAL(){
        @Override
        List<Double> getModifications(double p, double precursorMass) {
            double r = precursorMass - p - Consts.WATER;
            return Arrays.asList(p-1, p+1, p + Consts.CO, p + Consts.WATER, r-1, r+1, r - Consts.AMMONIA, r - Consts.WATER, r + Consts.WATER);
        }
    };

    abstract List<Double> getModifications(double p, double precursorMass);

}
