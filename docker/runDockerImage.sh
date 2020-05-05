docker run \
       -it \
       -v databases:/usr/local/databases \
       --mount type=bind,source="$(pwd)"/primerEvalOutput,target=/home/runDir/output \
       --name primerEval \
       matthewvc1/primer_eval:0.1.2 \
       bash
