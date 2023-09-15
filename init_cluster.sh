# CHANGE THIS DOT FILE DIRECTORY
DOT_FILE_DIR=$HOME/dotfile

# Follow jdhao install zsh and oh-my-zsh
#
# Installation zsh requires ncurses
wget ftp://ftp.gnu.org/gnu/ncurses/ncurses-6.1.tar.gz
tar xf ncurses-6.1.tar.gz
cd ncurses-6.1
./configure --prefix=$HOME/local CXXFLAGS="-fPIC" CFLAGS="-fPIC"
make -j && make install

# Build & install zsh
ZSH_SRC_NAME=$HOME/packages/zsh.tar.xz
ZSH_PACK_DIR=$HOME/packages/zsh
ZSH_LINK="https://sourceforge.net/projects/zsh/files/latest/download"

if [[ ! -d "$ZSH_PACK_DIR" ]]; then
    echo "Creating zsh directory under packages directory"
    mkdir -p "$ZSH_PACK_DIR"
fi

if [[ ! -f $ZSH_SRC_NAME ]]; then
    curl -Lo "$ZSH_SRC_NAME" "$ZSH_LINK"
fi

tar xJvf "$ZSH_SRC_NAME" -C "$ZSH_PACK_DIR" --strip-components 1
cd "$ZSH_PACK_DIR"

./configure --prefix="$HOME/local" \
    CPPFLAGS="-I$HOME/local/include" \
    LDFLAGS="-L$HOME/local/lib"
make -j && make install


echo <<EOF > ~/.bash_profile
export PATH=$HOME/local/bin:$PATH
export SHELL=`which zsh`
[ -f "$SHELL" ] && exec "$SHELL" -l
EOF

source ~/.bash_profile

# Follow gogoliri/misc (After install zsh & oh-my-zsh)
git clone https://github.com/romkatv/powerlevel10k $ZSH_CUSTOM/themes/powerlevel10k
git clone https://github.com/zsh-users/zsh-autosuggestions $ZSH_CUSTOM/plugins/zsh-autosuggestions
git clone https://github.com/zsh-users/zsh-syntax-highlighting $ZSH_CUSTOM/plugins/zsh-syntax-highlighting

# Backup zsh file and clone the zsh
cp ~/.zshrc ~/.zshrc_bku
cat $DOT_FILE_DIR/my_zshrc.sh > ~/.zshrc


# Get tmux source file
cp $DOT_FILE_DIR/my_tmux.conf ~/.tmux.conf
tmux source-file ~/.tmux.conf

# Install SpaceVim
curl -sLf https://spacevim.org/install.sh | bash


# Change default editor to vi (SpaceVim overwrite)
echo "export EDITOR=vi" >> ~/.zshrc


# Install Rust and Cargo
curl https://sh.rustup.rs -sSf | sh\n
source "$HOME/.cargo/env"

# Install ripgrep
cargo install ripgrep



# Install joshuto (ranger like)
cargo install --git https://github.com/kamiyaa/joshuto.git --force
echo "alias js=joshuto" >> ~/.zshrc

# Install spack [just gather NOT TEST, filter this out by order]
git clone -c feature.manyFiles=true https://github.com/spack/spack.git
. spack/share/spack/setup-env.sh

spack install gcc target=aarch64
spack location -i llvm@12

cd spack/opt/spack/linux-fedora37-neoverse_n1/gcc-12.2.1/llvm-12.0.1-ehp5a4hmurlr5x5rt4itka7ofpwsdjl6/bin/
