package bioinfo.energy.potential.preparation.voronoi;

/**
 * Datatype containing information about state of VoronoiData object
 * INIT is set after initialisation
 * CLOSE is set after writing all values into the object
 * CALC is set after calculatig the decomposition via VoroPPWrap and writing them into the object
 * @author andreseitz
 *
 */
public enum VoroDataState {

	INIT,CLOSE,CALC
}
