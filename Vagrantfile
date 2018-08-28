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
        aws.instance_type = "r5.12xlarge"
        aws.keypair_name = "forge"
        aws.subnet_id = "subnet-1fc8de7a"
        aws.security_groups = ["sg-38c9a872"]  # allows 22, 80 and 443
        aws.associate_public_ip = true
        aws.block_device_mapping = [{
            'DeviceName' => "/dev/sdf",
            'VirtualName' => "gp2_1",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdg",
            'VirtualName' => "gp2_2",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdh",
            'VirtualName' => "gp2_3",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdi",
            'VirtualName' => "gp2_4",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdj",
            'VirtualName' => "gp2_5",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdk",
            'VirtualName' => "gp2_6",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdl",
            'VirtualName' => "gp2_7",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        },
        {
            'DeviceName' => "/dev/sdm",
            'VirtualName' => "gp2_8",
            'Ebs.VolumeSize' => 50,
            'Ebs.DeleteOnTermination' => true,
            'Ebs.VolumeType' => 'gp2'
        }]
        override.ssh.username = "ec2-user"
        override.ssh.private_key_path = "~/.aws/forge.pem"
        # Good bids:
        #              vCPU  GiB mem  us-east-1  us-east-2
        # r5.4xlarge     16      128       0.35       0.20
        # r5.12xlarge    48      384       1.00       0.60
        # r5.24xlarge    96      768       1.90       1.10
        aws.region_config REGION do |region|
            region.spot_instance = true
            region.spot_max_price = "1.00"
        end
    end

    config.vm.provision "shell", privileged: true, name: "install Linux packages", inline: <<-SHELL
        yum install -q -y aws-cli wget unzip tree sysstat mdadm
    SHELL

    config.vm.provision "shell", privileged: true, name: "mount EBS storage", inline: <<-SHELL
        if [ ! -d /work ] ; then
            echo "Listing blocks:"
            lsblk
            # https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/raid-config.html
            
            echo "do RAID0 build"
            if [ ! -e "/dev/xvdf" ] ; then
                if [ ! -e "/dev/sdf" ] ; then
                    echo "**ERROR** -- no EBS drive at either /dev/sdf or /dev/xvdf"
                    exit 1
                fi
                mdadm --create --verbose /dev/md0 \
                    --level=0 --name=MY_RAID \
                    --raid-devices=8 /dev/sdf /dev/sdg /dev/sdh /dev/sdi \
                                     /dev/sdj /dev/sdk /dev/sdl /dev/sdm
            else
                mdadm --create --verbose /dev/md0 \
                    --level=0 --name=MY_RAID \
                    --raid-devices=8 /dev/xvdf /dev/xvdg /dev/xvdh /dev/xvdi \
                                     /dev/xvdj /dev/xvdk /dev/xvdl /dev/xvdm
            fi

            echo "wait for RAID0 build"
            cat /proc/mdstat
            sleep 10
            cat /proc/mdstat
            sleep 10
            cat /proc/mdstat
            sudo mdadm --detail /dev/md0
            
            echo "mkfs RAID0"
            mkfs -q -t ext4 -L MY_RAID /dev/md0
            
            echo "ensure RAID0 is reassembled automatically on boot"
            sudo mdadm --detail --scan | sudo tee -a /etc/mdadm.conf
            
            echo "Create ramdisk image to preload the block device modules"
            sudo dracut -H -f /boot/initramfs-$(uname -r).img $(uname -r)
            
            echo "Mount RAID0"
            mkdir /work
            mount LABEL=MY_RAID ${DRIVE} /work/
            chmod a+w /work
        fi
    SHELL

    config.vm.provision "file", source: "~/.aws/forge.pem", destination: "~ec2-user/.ssh/id_rsa"
    config.vm.provision "file", source: "~/.aws/credentials", destination: "~ec2-user/.aws/credentials"
    config.vm.provision "file", source: "~/.aws/config", destination: "~ec2-user/.aws/config"
    config.vm.provision "file", source: "./forge-contain.py", destination: "/work/forge-contain.py"

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
            --interpreter pypy3 \
            -- \
            --method hybrid \
            --reference /container-ref/hs37d5.fa \
            --vars /container-var/hs37d5.1ksnp \
            --cache-from /container-cache \
            --window-size 100 \
            --nslots 10 \
            --nthreads-per 6 \
            --prune 10 \
            --query-batch 40000000 \
            --output /container-output/hybrid.txt \
            --temp /container-temp
        
        echo "==Vagrantfile== Uploading rankings"
        aws s3 sync --quiet /container-output/ s3://forge-langmead/results/
    SHELL
end
