FROM continuumio/miniconda3:latest AS builder

WORKDIR /app

COPY environment.yml environment.yml
RUN conda env create -f environment.yml --name npmine_web_app
RUN apt-get update && apt-get install -y curl

SHELL ["conda", "run", "-n", "npmine_web_app", "/bin/bash", "-c"]

COPY . .

FROM continuumio/miniconda3:latest AS production

WORKDIR /app

COPY environment.yml environment.yml
RUN conda env create -f environment.yml --name npmine_web_app_prod

SHELL ["/bin/bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate npmine_web_app_prod && exec \"$@\""]

COPY populate_database.py populate_database.py 

EXPOSE 5000

CMD gunicorn --bind 0.0.0.0:5000 'websiteNPMINE:create_app()'
