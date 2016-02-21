(ns cuckoo-search-clj.test-blas-lapack
  (:import (no.uib.cipr.matrix DenseVector DenseMatrix SVD)))

(def cols (into-array DenseVector [(DenseVector. (double-array [1 3 3 3]))
                                   (DenseVector. (double-array [2 1 3 3]))
                                   (DenseVector. (double-array [2 2 1 3]))
                                   (DenseVector. (double-array [2 2 2 1]))]))

(def matA (DenseMatrix. cols))

(def svd (SVD. (.numRows matA) (.numColumns matA)))

(def svd' (.factor svd matA))

