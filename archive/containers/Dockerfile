FROM mambaorg/micromamba:1.5.8
COPY --chown=$MAMBA_USER:$MAMBA_USER environments/environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -a -y
WORKDIR /workspace/MembraneForger
COPY --chown=$MAMBA_USER:$MAMBA_USER . /workspace/MembraneForger
ENV PYTHONPATH=/workspace/MembraneForger
ENTRYPOINT ["./run_pipeline.sh"]
CMD ["doctor"]
