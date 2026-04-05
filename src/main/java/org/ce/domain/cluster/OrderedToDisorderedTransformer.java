package org.ce.domain.cluster;

import static org.ce.domain.cluster.ClusterPrimitives.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Transforms ordered-phase cluster coordinates to the disordered (HSP) reference frame.
 *
 * Applies a 3x3 rotation matrix and translation vector to all site positions
 * in the given clusters.
 */
public class OrderedToDisorderedTransformer {

    /**
     * Transforms a list of clusters using the given rotation matrix and translation vector.
     *
     * @param rotateMat 3x3 rotation matrix
     * @param translateMat 3-element translation vector
     * @param clusCoordCart clusters to transform
     * @return transformed clusters with new coordinates
     */
    public static List<Cluster> transform(
            double[][] rotateMat,
            double[] translateMat,
            List<Cluster> clusCoordCart) {

        List<Cluster> result = new ArrayList<>();

        for (Cluster cluster : clusCoordCart) {

            List<Sublattice> newSublattices = new ArrayList<>();

            for (Sublattice sub : cluster.getSublattices()) {

                List<Site> newSites = new ArrayList<>();

                for (Site site : sub.getSites()) {

                    Position r = site.getPosition();

                    double x =
                            rotateMat[0][0] * r.getX() +
                            rotateMat[0][1] * r.getY() +
                            rotateMat[0][2] * r.getZ() +
                            translateMat[0];

                    double y =
                            rotateMat[1][0] * r.getX() +
                            rotateMat[1][1] * r.getY() +
                            rotateMat[1][2] * r.getZ() +
                            translateMat[1];

                    double z =
                            rotateMat[2][0] * r.getX() +
                            rotateMat[2][1] * r.getY() +
                            rotateMat[2][2] * r.getZ() +
                            translateMat[2];

                    newSites.add(
                            new Site(
                                    new Position(x, y, z),
                                    site.getSymbol()
                            )
                    );
                }

                newSublattices.add(new Sublattice(newSites));
            }

            result.add(new Cluster(newSublattices));
        }

        return result;
    }
}




