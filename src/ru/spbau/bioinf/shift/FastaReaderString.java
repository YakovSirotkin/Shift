package ru.spbau.bioinf.shift;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;

public class FastaReaderString  {

    private ArrayList<Protein> proteins = new ArrayList<Protein>();

    public ArrayList<Protein> getProteins() {
        return proteins;
    }

    public FastaReaderString(File proteinDatabase)
            throws Exception {

        BufferedReader input = new BufferedReader(new InputStreamReader(new FileInputStream(
                proteinDatabase), "UTF-8"));

        String s;
        String cur = "";
        int proteinId = 0;
        while ((s = input.readLine()) != null) {
            if (s.startsWith(">")) {
                if (cur.length() > 0) {
                    proteins.add(new Protein(proteinId, cur));
                    proteinId++;
                    cur = "";
                }
            } else {
                cur += s.trim();
            }
        }
        proteins.add(new Protein(proteinId, cur));
    }
}
