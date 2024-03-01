FROM mambaorg/micromamba as conda

ARG CONDA_HOME=/opt/conda
ARG ENVIRON_NAME=workflow
ARG WF_HOME=/wf

RUN micromamba create -y \
    -c conda-forge -c bioconda \
    -n ${ENVIRON_NAME} \
    python "snakemake>=8.4.0" \
    numpy pandas biopython \
    click openpyxl
RUN ${CONDA_HOME}/envs/${ENVIRON_NAME}/bin/python -m pip install \
    "miniwdl>=1.11.0" \
    miniwdl-backend-bare


FROM memesuite/memesuite:5.5.5

ARG CONDA_HOME=/opt/conda
ARG ENVIRON_NAME=workflow
ARG WF_HOME=/wf

COPY . ${WF_HOME}
COPY --from=conda ${CONDA_HOME}/envs/${ENVIRON_NAME} ${CONDA_HOME}/envs/${ENVIRON_NAME}
ENV PATH="${CONDA_HOME}/envs/${ENVIRON_NAME}/bin:${WF_HOME}/workflow/scripts:$PATH"
ENV MINIWDL__SCHEDULER__CONTAINER_BACKEND=bare
