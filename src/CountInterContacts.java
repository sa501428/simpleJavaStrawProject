import javastraw.reader.Dataset;
import javastraw.reader.Matrix;
import javastraw.reader.basics.Chromosome;
import javastraw.reader.basics.ChromosomeHandler;
import javastraw.reader.block.Block;
import javastraw.reader.block.ContactRecord;
import javastraw.reader.mzd.MatrixZoomData;
import javastraw.reader.type.HiCZoom;
import javastraw.reader.type.NormalizationType;
import javastraw.tools.HiCFileTools;

import java.util.Arrays;
import java.util.List;

public class CountInterContacts {

    public static void main(String[] args){
        System.out.println("Args: "+Arrays.toString(args));
        if (args.length < 1){
            printUsageAndExit();
        }

        // Load the .hic file into a dataset
        Dataset ds = HiCFileTools.extractDatasetForCLT(args[0], false, false);

        // We can take a look at the resolutions available
        // List<HiCZoom> zooms = ds.getAllPossibleResolutions();
        // Or if you have a preferred resolution
        int resolution = 500000;
        HiCZoom zoom = ds.getZoomForBPResolution(resolution);

        // Let's pick our normalization
        // Typically you'd want SCALE or KR
        // e.g. NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString("SCALE");
        // (or) norm = ds.getNormalizationHandler().getNormTypeFromString("KR");
        // but for this example, I actually want no normalization (just counting)
        NormalizationType norm = ds.getNormalizationHandler().getNormTypeFromString("NONE");

        // Let's get the chromosomes
        ChromosomeHandler handler = ds.getChromosomeHandler();
        Chromosome[] chromosomes = handler.getChromosomeArrayWithoutAllByAll();

        // In this example, we are counting contacts
        long totalNumberOfIntraContacts = 0L;
        long totalNumberOfInterContacts = 0L;

        // iterating on the whole genome
        for(int i = 0; i < chromosomes.length; i++){
            for(int j = i; j < chromosomes.length; j++){

                // get the Matrix parent object
                Matrix matrix = ds.getMatrix(chromosomes[i], chromosomes[j]);
                if (matrix == null) continue;

                // let's get the MatrixZoomData object for a particular resolution
                MatrixZoomData zd = matrix.getZoomData(zoom);
                if (zd == null) continue;

                // now let's get the list of contact for a specific region
                // for this problem, we actually want to grab the whole chromosome region
                // so our bounding box limits will look like this:
                long rowStart = 0;
                long columnStart = 0;
                long rowEnd = chromosomes[i].getLength()/resolution;
                long columnEnd = chromosomes[j].getLength()/resolution;

                // Sometimes for Intra chromosomal regions you may want
                // to not just get the upper triangular matrix, but also
                // fill in the contacts below the diagonal
                // in such cases, fillUnderDiagonal should be set to true
                // but this isn't something we should do for inter regions
                boolean fillUnderDiagonal = false;
                List<Block> blocks = zd.getNormalizedBlocksOverlapping(rowStart, columnStart, rowEnd, columnEnd, norm,
                        false, fillUnderDiagonal);
                // iterate over the non-zero entries of the matrix
                long localSum = 0L;
                for(Block block : blocks){
                    for (ContactRecord cr : block.getContactRecords()){
                        // int row = cr.getBinX();
                        // int col = cr.getBinY();
                        float counts = cr.getCounts();
                        if(Float.isNaN(counts) || Float.isInfinite(counts)) continue;
                        localSum += counts;
                    }
                }
                if(i == j){
                    totalNumberOfIntraContacts += localSum;
                } else {
                    totalNumberOfInterContacts += localSum;
                }
            }
        }
        // correct for counts below diagonal
        totalNumberOfIntraContacts *= 2;
        totalNumberOfInterContacts *= 2;
        System.out.println("Total number of Intra-chromosomal contacts is "+totalNumberOfIntraContacts);
        System.out.println("Total number of Inter-chromosomal contacts is "+totalNumberOfInterContacts);
    }

    private static void printUsageAndExit() {
        System.err.println("Usage: java -jar counter.jar <inter.hic>");
        System.exit(2);
    }
}
