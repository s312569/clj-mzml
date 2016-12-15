(ns clj-mzml.core
  (:require [clojure.java.io :as io]
            [clojure.data.xml :refer [parse]]
            [clojure.zip :refer [node xml-zip]]
            [clojure.data.zip.xml :refer [xml-> attr xml1-> text attr]]
            [biodb.core :as bdb])
  (:import [org.apache.commons.codec.binary Base64]))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; utilities
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- parse-mzml
  [r]
  (let [x (parse r)]
    (if (= (:tag x) :indexedmzML)
      (first (:content x))
      x)))

(defn- get-parameter-tag
  [xml tag]
  (->> (take-while #(not (= (:tag %) :run))
                   (:content xml))
       (filter #(= (:tag %) tag))
       first
       xml-zip))

(defn- parse-cvs
  [z]
  (if (seq z)
    (->> (filter (fn [[k v]] (seq v))
                 {:cvparams (->> (map #(let [a (:attrs %)]
                                         (vector (:accession a) a))
                                      (xml-> z :cvParam node))
                                 (into {}))
                  :userparams (map :attrs (xml-> z :userParam node))
                  :references (->> (xml-> z :referenceableParamGroupRef (attr :ref)))})
         (into {}))))

(defn- parse-reference-group
  [zs]
  (->> (map #(vector (xml1-> % (attr :id))
                     (merge (-> (node %) :attrs (dissoc :id)) (parse-cvs %)))
            zs)
       (into {})))

(defn- components
  [z tag]
  (->> (xml-> z :componentList tag)
       (map #(merge {:order (xml1-> % (attr :order))} (parse-cvs %)))))

(defn- attr-cvs
  [z]
  (merge (-> z node :attrs) (parse-cvs z)))

(defn- inflate 
  [s]
  (let [h (new java.util.zip.InflaterInputStream
               (new java.io.ByteArrayInputStream s)
               (new java.util.zip.Inflater))]
    (try
      (loop [st h
             acc nil]
        (if (= (.available h) 0)
          (byte-array (map unchecked-byte (reverse acc)))
          (recur h (cons (.read h) acc))))
      (catch Exception e
        (println e))
      (finally
        (.close h)))))

(defn decode-binary
  [b]
  (let [decoded (.decode (Base64.) (:binary b))
        inflated (-> (if (get-in b [:cvparams "MS:1000576"]) decoded (inflate decoded))
                     java.nio.ByteBuffer/wrap
                     (.order java.nio.ByteOrder/LITTLE_ENDIAN))
        buffer (condp (fn [x y] (contains? y x)) (:cvparams b)
                 "MS:1000519" (.asIntBuffer inflated)
                 "MS:1000521" (.asFloatBuffer inflated)
                 "MS:1000523" (.asDoubleBuffer inflated)
                 "MS:1000522" (.asLongBuffer inflated)
                 :else (throw (Throwable. "Data type not supported")))]
    (loop [d buffer
           acc nil]
      (if (not (.hasRemaining d))
        (vec (reverse acc))
        (recur d (cons (.get d) acc))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; spectra
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- scans
  [z]
  (let [s (xml1-> z :scanList)]
    {:scans (merge (attr-cvs s)
                   {:scans (->> (xml-> s :scan)
                                (map #(merge (attr-cvs %)
                                             {:windows (->> (xml-> s :scanWindowList :scanWindow)
                                                            (map attr-cvs))})))})}))

(defn- precursors
  [z]
  {:precursors (map #(merge (-> % node :attrs)
                            {:isolation-window (-> (xml1-> % :isolationWindow)
                                                   parse-cvs)
                             :selected-ions (->> (xml-> % :selectedIonList :selectedIon)
                                                 (map parse-cvs))
                             :activation (-> (xml1-> % :activation)
                                             parse-cvs)})
                    (xml-> z :precursorList :precursor))})

(defn- products
  [z]
  {:products (map parse-cvs (xml-> z :productList :product :isolationWindow))})

(defn- binaries
  [z]
  {:binaries (->> (map #(merge (attr-cvs %) {:binary (xml1-> % :binary text)})
                       (xml-> z :binaryDataArrayList :binaryDataArray))
                  (map #(-> (assoc % :array (decode-binary %))
                            (dissoc :binary))))})

(defn- spec
  [z]
  (merge (attr-cvs z)
         (scans z)
         (precursors z)
         (products z)
         (binaries z)))

(defn spectra-seq
  [reader]
  (let [xml (parse-mzml reader)]
    (->> (filter #(#{:run} (:tag %)) (:content xml)) first :content
         (filter #(#{:spectrumList} (:tag %))) first :content
         (filter #(#{:spectrum} (:tag %)))
         (map xml-zip)
         (map spec))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; getting parameters, cv lists etc
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def parse-elements
  {:cvlist (fn [x] (->> (xml-> (get-parameter-tag x :cvList) :cv node) (map :attrs)))
   
   :file-content (fn [x] (->> (xml1-> (get-parameter-tag x :fileDescription) :fileContent)
                             parse-cvs))

   :source-files (fn [x]
                   (->> (xml-> (get-parameter-tag x :fileDescription) :sourceFileList :sourceFile)
                        parse-reference-group))

   :contacts (fn [x] (->> (xml-> (get-parameter-tag x :fileDescription) :contact)
                         (map parse-cvs)))

   :referenceable (fn [x] (->> (xml-> (get-parameter-tag x :referenceableParamGroupList)
                                     :referenceableParamGroup)
                              parse-reference-group))

   :samples (fn [x] (->> (xml-> (get-parameter-tag x :sampleList) :sample)
                        parse-reference-group))

   :software (fn [x] (->> (xml-> (get-parameter-tag x :softwareList) :software)
                         parse-reference-group))

   :scan-settings (fn [x] (->> (xml-> (get-parameter-tag x :scanSettingsList) :scan-settings)
                              (map #(vector (xml1-> % (attr :id))
                                            (merge (-> (node %) :attrs (dissoc :id))
                                                   (parse-cvs %)
                                                   {:source-refs (xml-> % :sourceFileRefList
                                                                        :sourceFileRef
                                                                        (attr :ref))
                                                    :targets (->> (xml-> % :targetList :target)
                                                                  (map parse-cvs))})))))
   :instrument-config (fn [x] (->> (xml-> (get-parameter-tag x :instrumentConfigurationList)
                                         :instrumentConfiguration)
                                  (map #(vector (xml1-> % (attr :id))
                                                (merge (-> (node %) :attrs (dissoc :id))
                                                       (parse-cvs %)
                                                       {:sources (components % :source)
                                                        :analyzers (components % :analyzer)
                                                        :detectors (components % :detector)})))
                                  (into {})))

   :data-processing (fn [x] (->> (xml-> (get-parameter-tag x :dataProcessingList) :dataProcessing)
                                (map #(vector (xml1-> % (attr :id))
                                              (map (fn [x]
                                                     (merge (-> x node :attrs)
                                                            (parse-cvs x)))
                                                   (xml-> % :processingMethod))))
                                (into {})))

   :default-data-proc-ref (fn [x] (->> (filter #(#{:run} (:tag %)) (:content x))
                                      first :content first :attrs :defaultDataProcessingRef))

   :spectrum-count (fn [x] (Integer/parseInt (->> (filter #(#{:run} (:tag %)) (:content x))
                                                 first :content first :attrs :count)))
   
   :run-parameter (fn [x]
                    (let [re (->> (filter #(#{:run} (:tag %)) (:content x)) first)
                          z (-> (assoc re :content
                                       (take-while #(not (#{:spectrumList :chromatogramList}
                                                          (:tag %)))
                                                   (:content re)))
                                xml-zip)]
                      (merge (-> z node :attrs) (parse-cvs z))))})

(defn mzml-parameters
  ([file] (mzml-parameters file false))
  ([file ks]
   (with-open [r (io/reader file)]
     (let [xml (parse-mzml r)]
       (->> (map (fn [[k v]] (vector k (v xml))) (if ks
                                                  (select-keys parse-elements ks)
                                                  parse-elements))
            (into {}))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; biodb compatability
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod bdb/table-spec :mzml
  [q]
  (vector [:accession :text "PRIMARY KEY"]
          [:number :integer "UNIQUE"]
          [:src :binary "NOT NULL"]))

(defmethod bdb/prep-sequences :mzml
  [q]
  (->> (:coll q)
       (map #(hash-map :accession (:id %) :number (Integer/parseInt (:index %))
                       :src (bdb/freeze %)))))

(defmethod bdb/restore-sequence :mzml
  [q]
  (bdb/thaw (:src (dissoc q :type))))
