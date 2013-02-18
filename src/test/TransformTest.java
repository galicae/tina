package test;

import bioinfo.proteins.PDBEntry;
import bioinfo.superpos.Kabsch;
import bioinfo.superpos.Transformation;
import cern.colt.*;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;

public class TransformTest {
	
	public static void main(String[] args) {
		double[][] arg1 = {{1,1,1},{2,3,2},{4,7,3}};
		DenseDoubleMatrix2D obj1 = new DenseDoubleMatrix2D(arg1);//rot
		
		double[][] arg2 = {{5,5,5},{8,8,8},{2,2,2}};
		DenseDoubleMatrix2D obj2 = new DenseDoubleMatrix2D(arg2);//grŸn
		
		double[][] arg3 = {{4,5,6},{3,2,1},{1,2,1}};
		DenseDoubleMatrix2D obj3 = new DenseDoubleMatrix2D(arg3);//schwarz
		
		System.out.println("obj1\n"+obj1);
		System.out.println("obj2\n"+obj2);
		System.out.println("obj3\n"+obj3);
		
		double[][][] pair21 = {arg1,arg2};
		Transformation trans21 = Kabsch.calculateTransformation(pair21);//grŸn -> rot
		DenseDoubleMatrix2D rot21 = trans21.peekRotation();
		DenseDoubleMatrix1D trl21 = trans21.peekTranslation();
		System.out.println("rot21\n"+rot21);
		System.out.println("trl21\n"+trl21);
		
		double[][][] pair12 = {arg2,arg1};
		Transformation trans12 = Kabsch.calculateTransformation(pair12);//grŸn -> rot
		DenseDoubleMatrix2D rot12 = trans12.peekRotation();
		DenseDoubleMatrix1D trl12 = trans12.peekTranslation();
		System.out.println("rot12\n"+rot12);
		System.out.println("trl12\n"+trl12);
		
		
		double[][][] pair32 = {arg2,arg3};
		Transformation trans32 = Kabsch.calculateTransformation(pair32);//schwarz -> grŸn
		DenseDoubleMatrix2D rot32 = trans32.peekRotation();
		DenseDoubleMatrix1D trl32 = trans32.peekTranslation();
		System.out.println("rot32\n"+rot32);
		System.out.println("trl32\n"+trl32.toString());

		
		
		
		
		DenseDoubleMatrix2D arg2onarg1 = new DenseDoubleMatrix2D(trans21.transform(arg2));//grŸn -> rot
		System.out.println(arg2onarg1.toString());
		DenseDoubleMatrix2D uni21 = createUniMatrix(rot21, trl21);
		DenseDoubleMatrix2D uni32 = createUniMatrix(rot32, trl32);
		System.out.println("uni21\n"+uni21.toString());
		System.out.println("uni32\n"+uni32.toString());
		
		Algebra alg = new Algebra();
		DenseDoubleMatrix2D result = (DenseDoubleMatrix2D)alg.mult(DoubleFactory2D.dense.identity(4),
				computeTransitiveTransformation((DenseDoubleMatrix2D)DoubleFactory2D.dense.identity(4), 
						uni32, uni21));
		
		System.out.println("action results\n"+result);
		System.out.println("moved points\n"+movePoints(result, obj3).toString());
		
		
	}
	
	public static DenseDoubleMatrix2D movePoints(DenseDoubleMatrix2D unifourxfour, DenseDoubleMatrix2D points){
		Algebra alg = new Algebra();
		DenseDoubleMatrix2D result = new DenseDoubleMatrix2D(points.rows(), 4);
		double tmp = 0.0d;
		DenseDoubleMatrix2D pointsfourxn = new DenseDoubleMatrix2D(points.rows(), 4);
		for(int i = 0; i != points.rows(); i++){
			for(int j = 0; j != 3; j++){
				pointsfourxn.set(i,j,points.get(i, j));
			}
			pointsfourxn.set(i, 3, 1);
			for(int j = 0; j != 4; j++){
				tmp = 0.0d;
				for(int k = 0; k != 4; k++){
					tmp += pointsfourxn.get(i, j)*unifourxfour.get(k, j);
				}
				result.set(i,j,tmp);
			}
		}
		return result;
//		return (DenseDoubleMatrix2D)alg.mult(pointsfourxn, unifourxfour);	
	}

	public static DenseDoubleMatrix2D createUniMatrix(DenseDoubleMatrix2D rot, DenseDoubleMatrix1D trl){
		DenseDoubleMatrix2D res= new DenseDoubleMatrix2D(4, 4);
			for(int i = 0; i != 3; i++){
				for(int j = 0; j != 3; j++){
					res.set(i, j, rot.get(i, j));
				}
			}
			for(int i = 0; i != 3; i++){
				res.set(i,3,trl.get(i));
			}
			for(int i = 0; i != 3; i++){
				res.set(3, i, 0);
			}
			res.set(3,3,1);
		return res;
	}
	
	public static DenseDoubleMatrix2D computeTransitiveTransformation(DenseDoubleMatrix2D priorTransformationOfMover, DenseDoubleMatrix2D actualMovementKabsch, DenseDoubleMatrix2D priorMovementOfBlueprint){
		Algebra alg = new Algebra();
		priorTransformationOfMover = (DenseDoubleMatrix2D)alg.inverse(priorTransformationOfMover);
		return (DenseDoubleMatrix2D)alg.mult(priorMovementOfBlueprint,alg.mult(priorTransformationOfMover, actualMovementKabsch));		
//		return (DenseDoubleMatrix2D)alg.mult(alg.mult(actualMovementKabsch,priorTransformationOfMover),priorMovementOfBlueprint);		

	
	}
}


