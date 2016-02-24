;Author Rigoberto Leander Salgado Reyes <rlsalgado2006 @gmail.com>
;
;Copyright 2016 by Rigoberto Leander Salgado Reyes.
;
;This program is licensed to you under the terms of version 3 of the
;GNU Affero General Public License. This program is distributed WITHOUT
;ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
;MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
;AGPL (http:www.gnu.org/licenses/agpl-3.0.txt) for more details.

(ns cuckoo-search-clj.functions
  ;  (:import (no.uib.cipr.matrix DenseMatrix SVD))
  (:import (org.netlib.lapack Dgesvd))
  (:import (org.netlib.util intW)))

(defn make-column
  ""
  [seed size]
  (concat (repeat (- (count seed) size) 0) (subvec seed 0 size)))

(defn make-toeplitz-matrix
  "Toeplitz Matrix in column major form"
  [seed]
  (double-array (mapcat (partial make-column seed) (range (count seed) 0 -1))))

;(defn svd
;  ""
;  [seed matrix-maker vectors?]
;  (let [size (count seed)]
;    (.factor (SVD. size size vectors?) (DenseMatrix. size size (matrix-maker seed) false))))

(defn svd
  ""
  [seed matrix-maker vectors?]
  (let [size (count seed)
        s (double-array size 0.0)
        u (double-array (* size size) 0.0)
        vt (double-array (* size size) 0.0)
        work (double-array (* size size 2) 0.0)
        info (intW. 2)
        get-vec (if vectors? "A" "N")
        _ (Dgesvd/dgesvd get-vec get-vec size size (matrix-maker seed) 0 size s 0 u 0 size vt 0 size
            work 0 (* size size 2) info)]
    {:s s :u u :vt vt :info info}))

(def calc-sv #(:s (svd % %2 false)))

(def zeros #(vec (repeat (* % %2) 0.0)))

(defn set-diag
  ""
  ([matrix fun rows row]
    (set-diag matrix fun rows row rows 0))
  ([matrix fun rows row cols col]
    (if (>= (max row col) (min rows cols)) matrix
      (recur (update matrix (+ row (* col cols)) fun) fun rows (inc row) cols (inc col)))))

(defn ddot
  "Dot product"
  ([coll1 coll2]
    (ddot coll1 coll2 0))
  ([coll1 coll2 acc]
    (if (or (empty? coll1) (empty? coll2)) acc
      (recur (rest coll1) (rest coll2) (+ acc (* (first coll1) (first coll2)))))))

(defn dgemv
  "Matrix-Vector multiplication"
  ([matrix rows cols vector]
    (dgemv matrix rows cols vector (vec (repeat rows 0.0)) 0))
  ([matrix rows cols vector result pos]
    (if (or (>= pos (* rows cols)) (not= (count matrix) (* rows cols)) (not= (count vector) cols))
      result
      (recur matrix rows cols vector
        (update result (rem pos rows) #(+ % (* (get matrix pos) (get vector (quot pos rows)))))
        (inc pos)))))

(defn fill-J
  ""
  [J Ai qi pi row col size]
  (if (>= col size) J
    (recur
      (assoc J (+ row (* col size))
        (ddot pi (dgemv (set-diag Ai inc size col) size size qi)))
      Ai qi pi row (inc col) size)))

(defn JacIASVPToeplitzTriInf
  ""
  ([seed matrix-maker]
    (let [{P :u Q :vt} (svd seed matrix-maker true)
          size (count seed)
          Ai (zeros size size)
          J (zeros size size)]
      (JacIASVPToeplitzTriInf (vec P) (vec Q) Ai J 0 size)))
  ([P Q Ai J row size]
    (let [init (* row size)
          get-col #(subvec % init (+ init size))]
      (if (>= row size) J
        (recur P Q Ai (fill-J J Ai (get-col Q) (get-col P) row 0 size) (inc row) size)))))

(def norm2 #(reduce (fn [acc v] (+ acc (* v v))) 0.0 %))
;
;(defn descending
;  ""
;  [n2fnewx n2fx newx stop-tol seed fun]
;  (if (<= (- n2fnewx n2fx) stop-tol) newx
;    (let [newx (map #(/ (+ % %2) 2) newx seed)]
;      (recur (->> newx fun norm2) n2fx newx stop-tol seed fun))))
;
;(defn newtonBiseccionNLES
;  ""
;  [fun seed Jac tol abs-tol max-it]
;  1)

























