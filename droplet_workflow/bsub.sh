mkdir -p bsub_logs
bsub -q ser -n 2 -M 4000 -R 'span[hosts=1]' -o bsub_logs/run.log bash submit_snake.sh
