;Author Rigoberto Leander Salgado Reyes <rlsalgado2006 @gmail.com>
;
;Copyright 2016 by Rigoberto Leander Salgado Reyes.
;
;This program is licensed to you under the terms of version 3 of the
;GNU Affero General Public License. This program is distributed WITHOUT
;ANY EXPRESS OR IMPLIED WARRANTY, INCLUDING THOSE OF NON-INFRINGEMENT,
;MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. Please refer to the
;AGPL (http:www.gnu.org/licenses/agpl-3.0.txt) for more details.

(ns cuckoo-search-clj.operators)

(def pow #(Math/pow % %2))
(def sin #(Math/sin %))
(def abs #(Math/abs %))
(def exp #(Math/exp %))
(def PI Math/PI)
(def isNaN #(Double/isNaN %))

(declare gamma check-bounds get-cuckoos-process-problem check-best-nest better-fitness empty-nest-process-problem)

;------------------------------------------------------------------------------------------------
;GET-CUCKOOS
;------------------------------------------------------------------------------------------------
(defn get-cuckoos
  "Get cuckoos by ramdom walk"
  [{nest :nest, :as cs}]
  (let [beta (/ 3.0 2.0)
        sigma (pow (/ (* (gamma (+ 1.0 beta)) (sin (/ (* PI beta) 2.0)))
                     (* (gamma (/ (+ 1.0 beta) 2.0)) beta (pow 2.0 (/ (- beta 1.0) 2.0)))) (/ 1.0 beta))]
    (assoc cs :new-nest (mapv (partial get-cuckoos-process-problem cs beta sigma) nest))))

(defn get-cuckoos-process-solution
  ""
  [cs {best-solution :solution} beta sigma solution]
  (let [u-j (* (rand) sigma)
        v-j (rand)
        step-j (pow (/ u-j (abs v-j)) (/ 1.0 beta))]
    (mapv #(check-bounds cs (+ % (* (* 0.01 step-j (- % %2)) (rand)))) solution best-solution)))

(defn get-cuckoos-process-problem
  ""
  [{:keys [best-nest fun], :as cs} beta sigma p]
  (->> p
    (#(update % :solution (partial get-cuckoos-process-solution cs best-nest beta sigma)))
    (#(assoc % :fitness (fun %)))))

(defn check-bounds
  ""
  [{:keys [lb ub]} value]
  (if (isNaN value)
    lb
    (max (min value ub) lb)))

(defn gamma
  "implements Weirstrass's form (infinite product)
   Mathematical methods for Physicists, 4th ed. page 594."
  ([z] (gamma z 1000))
  ([z nterm] (let [g 0.577216 ; Euler-Mascheroni constant
                   ret-val (* z (exp (* g z)))
                   fn #(* % (* (+ 1 (/ z %2)) (exp (/ (- z) %2))))]
               (/ 1.0 (reduce fn ret-val (range 1 (inc nterm)))))))


;------------------------------------------------------------------------------------------------
;BEST-NEST
;------------------------------------------------------------------------------------------------
(defn best-nest
  "Find the current best nest"
  [{:keys [nest new-nest], :as cs}]
  (->> cs
    (#(assoc % :nest (mapv better-fitness nest new-nest)))
    check-best-nest))

(defn check-best-nest
  ""
  [{:keys [best-nest nest], :as cs}]
  (assoc cs :best-nest (reduce better-fitness best-nest nest)))

(defn better-fitness
  ""
  [{fitness1 :fitness, :as problem1} {fitness2 :fitness, :as problem2}]
  (if (<= fitness2 fitness1) problem2 problem1))


;------------------------------------------------------------------------------------------------
;EMPTY-NEST
;------------------------------------------------------------------------------------------------
(defn empty-nest
  "Replace some nests by constructing new solutions/nests"
  [cs]
  (let [rand* (rand)
        {:keys [nest perm1 perm2], :as cs} (->> cs
                                             (#(update % :perm1 shuffle))
                                             (#(update % :perm2 shuffle)))]
    (assoc cs :new-nest (mapv (partial empty-nest-process-problem rand* cs) nest perm1 perm2))))

(defn empty-nest-process-solution
  [rand* {nest :nest} pos1 pos2 solution]
  (let [solution-perm1 (get-in nest [pos1 :solution])
        solution-perm2 (get-in nest [pos2 :solution])]
    (mapv #(+ % (* rand* (- %2 %3))) solution solution-perm1 solution-perm2)))

(defn empty-nest-process-problem
  ""
  [rand* {:keys [pa fun], :as cs} problem pos-perm1 pos-perm2]
  (if (> (rand) pa)
    (->> problem
      (#(update % :solution (partial empty-nest-process-solution rand* cs pos-perm1 pos-perm2)))
      (#(assoc % :fitness (fun %))))
    problem))
;------------------------------------------------------------------------------------------------
