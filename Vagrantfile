# -*- mode: ruby -*-
# vi: set ft=ruby :

# Steps:
# 1. (install vagrant)
# 2. vagrant plugin install vagrant-aws-mkubenka --plugin-version "0.7.2.pre.22"
# 3. vagrant box add dummy https://github.com/mitchellh/vagrant-aws/raw/master/dummy.box
#
# Note: the standard vagrant-aws plugin does not have spot support

ENV['VAGRANT_DEFAULT_PROVIDER'] = 'aws'
REGION = "us-east-1"

Vagrant.configure("2") do |config|

    config.vm.box = "dummy"
    config.vm.synced_folder ".", "/vagrant", disabled: true

    config.vm.provider :aws do |aws, override|
        aws.region = REGION
        aws.ami = "ami-13401669"
        aws.tags = { 'Application' => 'forge' }
        aws.instance_type = "c5.4xlarge"
        aws.keypair_name = "forge"
        aws.subnet_id = "subnet-1fc8de7a"
        aws.security_groups = ["sg-38c9a872"]  # allows 22, 80 and 443
        aws.associate_public_ip = true
        aws.block_device_mapping = [{
            'DeviceName' => "/dev/sdf",
            'VirtualName' => "ephemeral0",
            'Ebs.VolumeSize' => 100,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        }]
        override.ssh.username = "ec2-user"
        override.ssh.private_key_path = "~/.aws/forge.pem"
        aws.region_config REGION do |region|
            region.spot_instance = true
            region.spot_max_price = "0.40"
        end
    end

    config.vm.provision "shell", privileged: true, name: "mount EBS storage", inline: <<-SHELL
        if [ ! -d /work ] ; then
            DRIVE=/dev/xvdf
            if [ ! -e "${DRIVE}" ] ; then
                DRIVE=/dev/sdf
            fi
            if [ ! -e "${DRIVE}" ] ; then
                echo "**ERROR** -- no EBS drive at either /dev/sdf or /dev/xvdf"
                exit 1
            fi
            mkfs -q -t ext4 ${DRIVE}
            mkdir /work
            mount ${DRIVE} /work/
            chmod a+w /work
        fi
    SHELL

    config.vm.provision "file", source: "~/.aws/forge.pem", destination: "~ec2-user/.ssh/id_rsa"
    config.vm.provision "file", source: "~/.aws/credentials", destination: "~ec2-user/.aws/credentials"
    config.vm.provision "file", source: "~/.aws/config", destination: "~ec2-user/.aws/config"
    config.vm.provision "file", source: "./forge-contain.py", destination: "/work/forge-contain.py"

    config.vm.provision "shell", privileged: true, name: "install Linux packages", inline: <<-SHELL
        yum install -q -y aws-cli wget unzip tree
    SHELL

    config.vm.provision "shell", privileged: false, name: "download inputs", inline: <<-SHELL
        mkdir -p /work/kmc_dbs
        echo "==Vagrantfile== Syncing phasing information"
        aws s3 sync --quiet s3://forge-langmead/input/phasing/ /work/phasing/
        echo "==Vagrantfile== Syncing references"
        aws s3 sync --quiet s3://forge-langmead/input/ref/     /work/ref/
        echo "==Vagrantfile== Syncing variants"
        aws s3 sync --quiet s3://forge-langmead/input/variant/ /work/variant/
        echo "==Vagrantfile== Syncing KMC DBs"
        aws s3 sync --quiet s3://forge-langmead/kmc_dbs/       /work/kmc_dbs/
        echo "Tree:"
        tree /work
    SHELL

    config.vm.provision "shell", privileged: false, name: "docker run forge", inline: <<-SHELL
        echo "==Vagrantfile== Running FORGe"
        mkdir -p /work/output
        sudo python /work/forge-contain.py \
            --mount /work/temp:/container-temp \
                    /work/kmc_dbs:/container-cache \
                    /work/ref:/container-ref \
                    /work/variant:/container-var \
                    /work/phasing:/container-phasing \
                    /work/output:/container-output \
            -- \
            --method hybrid \
            --reference /container-ref/hs37d5.fa \
            --vars /container-var/hs37d5.1ksnp \
            --cache-from /container-cache \
            --window-size 100 \
            --threads 16 \
            --prune 15 \
            --output /container-output/hybrid.txt \
            --temp /container-temp
        
        echo "==Vagrantfile== Uploading rankings"
        aws s3 sync --quiet /container-output/ s3://forge-langmead/results/
    SHELL
end
